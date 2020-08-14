/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016-2020 hyStrath
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of hyStrath, a derivative work of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "MillikanWhite.H"
#include "addToRunTimeSelectionTable.H"

#include "DynamicList.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace VTRelaxationModels
    {
        defineTypeNameAndDebug(MillikanWhite, 0);
        addToRunTimeSelectionTable
        (
            VTRelaxationModel,
            MillikanWhite,
            dictionary
        );
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::VTRelaxationModels::MillikanWhite::MillikanWhite
(
    const word& name1,
    const word& name2,
    const label& lname1,
    const label& lname2,
    const dictionary& dict2T,
    const dictionary& dictThermoPhy,
    const volScalarField& p,
    const volScalarField& Tt,
    const PtrList<volScalarField>& Tv,
    const PtrList<volScalarField>& nD
)
:
    VTRelaxationModel(name1, name2, lname1, lname2, dict2T, dictThermoPhy, p, Tt, Tv, nD)
{
    const scalar W1 = readScalar(dictThermoPhy.subDict(name1).subDict("specie").lookup("molWeight"));
    const scalar W2 = readScalar(dictThermoPhy.subDict(name2).subDict("specie").lookup("molWeight"));
    DynamicList<scalar> vibData(dictThermoPhy.subDict(name1).subDict("thermodynamics").lookup("vibrationalList"));
    const scalar TH1 = vibData[1];

    word subDictName = word::null;

    if (not VTFullCoeffsForm_)
    {
        scalar W12 = (W1 * W2) / (W1 + W2);
        A12_ = sqrt(W12) * pow(TH1, 4.0/3.0);
        B12_ = pow(W12, 0.25);
        scalar preAij = 0.0;
        scalar preMij = 0.0;

        if (not VTOverwriteDefault_)
        {
            preAij  = 1.16e-3;
            preMij  = 0.015;
            offset_ = 18.42;
        }
        else if (VTSpeciesDependent_ and VTCollidingPartner_)
        {
            if (dict2T.subDict("MillikanWhiteCoefficients").isDict(name1+"_"+name2))
            {
                subDictName = name1+"_"+name2;
            }
            else if (dict2T.subDict("MillikanWhiteCoefficients").isDict(name2+"_"+name1))
            {
                subDictName = name2+"_"+name1;
            }
            else if (dict2T.subDict("MillikanWhiteCoefficients").isDict(name1))
            {
                subDictName = name1;
            }
            else
            {
                subDictName = "allSpecies";

            }
        }
        else if (VTSpeciesDependent_ and dict2T.subDict("MillikanWhiteCoefficients").isDict(name1))
        {
            subDictName = name1;
        }
        else
        {
            subDictName = "allSpecies";
        }

        preAij = readScalar(dict2T.subDict("MillikanWhiteCoefficients").subDict(subDictName).lookup("preAij"));
        preMij = readScalar(dict2T.subDict("MillikanWhiteCoefficients").subDict(subDictName).lookup("preMij"));
        A12_ *= preAij;
        B12_ *= preMij;
    }
    else
    {
        if (not VTOverwriteDefault_)
        {
            A12_  = 221.53;
            B12_  = 0.029;
            offset_ = 18.42;
        }
        else if (VTSpeciesDependent_ and VTCollidingPartner_)
        {
            if (dict2T.subDict("MillikanWhiteCoefficients").isDict(name1+"_"+name2))
            {
                subDictName = name1+"_"+name2;
            }
            else if (dict2T.subDict("MillikanWhiteCoefficients").isDict(name2+"_"+name1))
            {
                subDictName = name2+"_"+name1;
            }
            else if (dict2T.subDict("MillikanWhiteCoefficients").isDict(name1))
            {
                subDictName = name1;
            }
            else
            {
                subDictName = "allSpecies";
            }
        }
        else if (VTSpeciesDependent_ and dict2T.subDict("MillikanWhiteCoefficients").isDict(name1))
        {
            subDictName = name1;
        }
        else
        {
            subDictName = "allSpecies";
        }
    }

    offset_ = readScalar(dict2T.subDict("MillikanWhiteCoefficients").subDict(subDictName).lookup("offset"));

}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::VTRelaxationModels::MillikanWhite::tauVT() const
{
    const fvMesh& mesh = this->T_.mesh();

    tmp<volScalarField> ttauVT
    (
        new volScalarField
        (
            IOobject
            (
                "tauVT_" + name1_ + "_" + name2_,
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimTime
        )
    );

    volScalarField& tauVT = ttauVT.ref();

    forAll(this->T_, celli)
    {
        tauVT[celli] =
            1.01325e5 / this->p_[celli] * exp(A12_*(pow(this->T_[celli], -1.0/3.0)
                - B12_) - offset_);
    }


    forAll(this->T_.boundaryField(), patchi)
    {
        const fvPatchScalarField& pTt = this->T_.boundaryField()[patchi];
        const fvPatchScalarField& pp = this->p_.boundaryField()[patchi];
        fvPatchScalarField& ptauVT = tauVT.boundaryFieldRef()[patchi];

        forAll(pTt, facei)
        {
            ptauVT[facei] =
            1.01325e5 / pp[facei] * exp(A12_*(pow(pTt[facei], -1.0/3.0)
                - B12_) - offset_);
        }
    }

    return ttauVT;
}


Foam::tmp<Foam::scalarField> Foam::VTRelaxationModels::MillikanWhite::tauVT
(
    const label patchi,
    const scalarField& p,
    const scalarField& Tt,
    const PtrList<scalarField>& Tv,
    const PtrList<scalarField>& nD
) const
{
    tmp<scalarField> ttauVT(new scalarField(Tt.size()));
    scalarField& tauVT = ttauVT.ref();

    forAll(Tt, facei)
    {
        tauVT[facei] =
            1.01325e5 / p[facei] * exp(A12_*(pow(Tt[facei], -1.0/3.0)
                - B12_) - offset_);
    }

    return ttauVT;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
