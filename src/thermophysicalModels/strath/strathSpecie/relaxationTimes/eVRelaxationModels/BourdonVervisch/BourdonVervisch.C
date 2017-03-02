/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

#include "BourdonVervisch.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace eVRelaxationModels
    {
        defineTypeNameAndDebug(BourdonVervisch, 0);
        addToRunTimeSelectionTable
        (
            eVRelaxationModel,
            BourdonVervisch, 
            dictionary
        );
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::eVRelaxationModels::BourdonVervisch::BourdonVervisch
(
    const word& name1,
    const label& lname1,
    const dictionary& dict2T,
    const dictionary& dictThermoPhy,
    const volScalarField& p,
    const volScalarField& Tv
)
:
    eVRelaxationModel(name1, lname1, dict2T, dictThermoPhy, p, Tv)
{   
    species1_ = lname1; name1_ = name1;
    
    Tlow_ = readScalar(dict2T.subDict("BourdonVervischCoefficients").lookup("Tlow"));
    Tmed_ = readScalar(dict2T.subDict("BourdonVervischCoefficients").lookup("Tmed"));
    Thigh_ = readScalar(dict2T.subDict("BourdonVervischCoefficients").lookup("Thigh"));
    
    word subDictName = word::null;
    
    if (not eVOverwriteDefault_)
    {
        BV1_[0]  = 5.019; BV1_[1] = 2.448;
        BV2_[0]  = -38.025; BV2_[1] = -18.704;
        BV3_[0]  = 64.219; BV3_[1] = 25.635;
    }
    else 
    {
        if (eVSpeciesDependent_ and dict2T.subDict("BourdonVervischCoefficients").isDict(name1))
        {        
            subDictName = name1; 
        }    
        else
        {
            subDictName = "allSpecies";
        }    
        
        BV1_ = dict2T.subDict("BourdonVervischCoefficients").subDict(subDictName).lookup("BV1");
        BV2_ = dict2T.subDict("BourdonVervischCoefficients").subDict(subDictName).lookup("BV2");
        BV3_ = dict2T.subDict("BourdonVervischCoefficients").subDict(subDictName).lookup("BV3");
    }
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::eVRelaxationModels::BourdonVervisch::taueV() const
{
    const fvMesh& mesh = this->Tv_.mesh();

    tmp<volScalarField> ttaueV
    (
        new volScalarField
        (
            IOobject
            (
                "taueV_" + name1_,
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimTime
        )
    );

    volScalarField& taueV = ttaueV();
    
    forAll(this->Tv_, celli)
    {
        //partial pressure of electrons greater than a given value = electrons in the cell
        if(this->p_[celli] < 1e-2)
        {
            taueV[celli] = Foam::GREAT;
        }
        else
        {
            FixedList<scalar, 2> tempRange = findTemperatureRange(this->Tv_[celli]);
            scalar BV1 = boundedLinearInterpolation(this->Tv_[celli], tempRange[0], tempRange[1], BV1_[0], BV1_[1]);
            scalar BV2 = boundedLinearInterpolation(this->Tv_[celli], tempRange[0], tempRange[1], BV2_[0], BV2_[1]);
            scalar BV3 = boundedLinearInterpolation(this->Tv_[celli], tempRange[0], tempRange[1], BV3_[0], BV3_[1]);
            
            //Info << log10(this->Tv_[celli]) << tab << BV1 << tab << BV2 << tab << BV3 << tab << BV1*sqr(log10(this->Tv_[celli])) + BV2*log10(this->Tv_[celli]) + BV3 << endl;
            taueV[celli] = 1.01325e5 / this->p_[celli] * pow(10, BV1*sqr(log10(this->Tv_[celli])) + BV2*log10(this->Tv_[celli]) + BV3);
        }
    }
    

    forAll(this->Tv_.boundaryField(), patchi)
    {
        const fvPatchScalarField& pTv = this->Tv_.boundaryField()[patchi];
        const fvPatchScalarField& pp = this->p_.boundaryField()[patchi];
        fvPatchScalarField& ptaueV = taueV.boundaryField()[patchi];

        forAll(pTv, facei)
        {
            if(pp[facei] < 1e-2)
            {
                ptaueV[facei] = Foam::GREAT;
            }
            else
            {
                FixedList<scalar, 2> tempRange = findTemperatureRange(pTv[facei]);
                scalar pBV1 = boundedLinearInterpolation(pTv[facei], tempRange[0], tempRange[1], BV1_[0], BV1_[1]);
                scalar pBV2 = boundedLinearInterpolation(pTv[facei], tempRange[0], tempRange[1], BV2_[0], BV2_[1]);
                scalar pBV3 = boundedLinearInterpolation(pTv[facei], tempRange[0], tempRange[1], BV3_[0], BV3_[1]);
                
                ptaueV[facei] = 1.01325e5 / pp[facei] * pow(10, pBV1*sqr(log10(pTv[facei])) + pBV2*log10(pTv[facei]) + pBV3);
            }
        }
    }

    return ttaueV;
}


Foam::tmp<Foam::scalarField> Foam::eVRelaxationModels::BourdonVervisch::taueV
(
    const label patchi,
    const scalarField& p,
    const scalarField& Tv
) const
{
    tmp<scalarField> ttaueV(new scalarField(Tv.size()));
    scalarField& taueV = ttaueV();
    
    forAll(Tv, facei)
    {
        if(p[facei] < 1e-2)
        {
            taueV[facei] = Foam::GREAT;
        }
        else
        {
            FixedList<scalar, 2> tempRange = findTemperatureRange(Tv[facei]);
            scalar pBV1 = boundedLinearInterpolation(Tv[facei], tempRange[0], tempRange[1], BV1_[0], BV1_[1]);
            scalar pBV2 = boundedLinearInterpolation(Tv[facei], tempRange[0], tempRange[1], BV2_[0], BV2_[1]);
            scalar pBV3 = boundedLinearInterpolation(Tv[facei], tempRange[0], tempRange[1], BV3_[0], BV3_[1]);
            
            taueV[facei] = 1.01325e5 / p[facei]* pow(10, pBV1*sqr(log10(Tv[facei])) + pBV2*log10(Tv[facei]) + pBV3);
        }
    }

    return ttaueV;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
