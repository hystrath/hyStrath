/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
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

#include "Knab.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace VVRelaxationModels
    {
        defineTypeNameAndDebug(Knab, 0);
        addToRunTimeSelectionTable
        (
            VVRelaxationModel,
            Knab, 
            dictionary
        );
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::VVRelaxationModels::Knab::Knab
(
    const word& name1,
    const word& name2,
    const label& lname1,
    const label& lname2,
    const dictionary& dict1,
    const dictionary& dict2,
    const volScalarField& p,
    const volScalarField& Tt,
    const PtrList<volScalarField>& Tv,
    const PtrList<volScalarField>& nD
)
:
    VVRelaxationModel(name1, name2, lname1, lname2, dict1, dict2, p, Tt, Tv, nD)
{   
    W1_ = 1.0e-3*readScalar(dict2.subDict(name1).subDict("specie").lookup("molWeight"));
    const scalar W2 = 1.0e-3*readScalar(dict2.subDict(name2).subDict("specie").lookup("molWeight"));
      
    W12_ = (W1_ * W2) / (W1_ + W2);
    
    if (not VVOverwriteDefault_)
    {
        P21_ = 0.01;
        sigma12_ = 1.0e-20;
    }
    else if (VVSpeciesDependent_ and VVCollidingPartner_)
    {        
        if (dict1.subDict("KnabCoefficients").isDict(name1+"_"+name2))
        {
            P21_ = readScalar(dict1.subDict("KnabCoefficients").subDict(name1+"_"+name2).lookup("P21"));
            sigma12_ = readScalar(dict1.subDict("KnabCoefficients").subDict(name1+"_"+name2).lookup("sigma12"));
        }
        else if (dict1.subDict("KnabCoefficients").isDict(name2+"_"+name1))
        {
            P21_ = readScalar(dict1.subDict("KnabCoefficients").subDict(name2+"_"+name1).lookup("P21"));
            sigma12_ = readScalar(dict1.subDict("KnabCoefficients").subDict(name2+"_"+name1).lookup("sigma12"));
        }
        else if (dict1.subDict("KnabCoefficients").isDict(name1))
        {
            P21_ = readScalar(dict1.subDict("KnabCoefficients").subDict(name1).lookup("P21"));
            sigma12_ = readScalar(dict1.subDict("KnabCoefficients").subDict(name1).lookup("sigma12"));
        }
        else
        {
            P21_ = readScalar(dict1.subDict("KnabCoefficients").subDict("allSpecies").lookup("P21"));
            sigma12_ = readScalar(dict1.subDict("KnabCoefficients").subDict("allSpecies").lookup("sigma12"));
        }    
    }
    else if (VVSpeciesDependent_ and dict1.subDict("KnabCoefficients").isDict(name1))
    {        
        P21_ = readScalar(dict1.subDict("KnabCoefficients").subDict(name1).lookup("P21"));
        sigma12_ = readScalar(dict1.subDict("KnabCoefficients").subDict(name1).lookup("sigma12"));
    } 
    else
    {
        P21_ = readScalar(dict1.subDict("KnabCoefficients").subDict("allSpecies").lookup("P21"));
        sigma12_ = readScalar(dict1.subDict("KnabCoefficients").subDict("allSpecies").lookup("sigma12"));   
    }
    
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::VVRelaxationModels::Knab::tauVV() const
{
    const fvMesh& mesh = this->Tt_.mesh();

    tmp<volScalarField> ttauVV
    (
        new volScalarField
        (
            IOobject
            (
                "tauVV_" + name1_ + "_" + name2_,
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionSet(0, 0, 1, 0, 0)
        )
    );

    volScalarField& tauVV = ttauVV.ref();

    forAll(this->Tt_, celli)
    {
        tauVV[celli] = P21_*sigma12_/sqrt(W12_);
    }
    

    forAll(this->Tt_.boundaryField(), patchi)
    {
        fvPatchScalarField& ptauVV = tauVV.boundaryFieldRef()[patchi];

        forAll(ptauVV, facei)
        {
            ptauVV[facei] = P21_*sigma12_/sqrt(W12_);
        }
    }

    return ttauVV;
}


Foam::tmp<Foam::scalarField> Foam::VVRelaxationModels::Knab::tauVV
(
    const label patchi,
    const scalarField& p,
    const scalarField& Tt,
    const PtrList<scalarField>& Tv,
    const PtrList<scalarField>& nD
) const
{
    tmp<scalarField> ttauVV(new scalarField(Tt.size()));
    scalarField& tauVV = ttauVV.ref();

    forAll(Tt, facei)
    {
        tauVV[facei] = P21_*sigma12_/sqrt(W12_);
    }

    return ttauVV;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
