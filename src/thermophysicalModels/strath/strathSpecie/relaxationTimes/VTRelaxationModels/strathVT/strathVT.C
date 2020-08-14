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

#include "strathVT.H"
#include "addToRunTimeSelectionTable.H"

#include "DynamicList.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace VTRelaxationModels
    {
        defineTypeNameAndDebug(strathVT, 0);
        addToRunTimeSelectionTable
        (
            VTRelaxationModel,
            strathVT,
            dictionary
        );
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::VTRelaxationModels::strathVT::strathVT
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
    species1_ = lname1;
    W1_ = readScalar(dictThermoPhy.subDict(name1).subDict("specie").lookup("molWeight"));
    const scalar W2 = readScalar(dictThermoPhy.subDict(name2).subDict("specie").lookup("molWeight"));
    DynamicList<scalar> vibData(dictThermoPhy.subDict(name1).subDict("thermodynamics").lookup("vibrationalList"));
    THETA1_ = vibData[1];

    if (not VTFullCoeffsForm_)
    {
        scalar W12 = (W1_ * W2) / (W1_ + W2);
        A12_ = sqrt(W12) * pow(THETA1_, 4.0/3.0);
        B12_ = pow(W12, 0.25);
        scalar preAij = 0.0;
        scalar preMij = 0.0;

        if (not VTOverwriteDefault_)
        {
            preAij  = 1.16e-3;
            preMij  = 0.015;
            offset_ = 18.42;
            sigma1_ = 1.0e-21;
            sigma2_ = 5.0e4;

            LambdaD_ = 0.0;
            LambdaE_ = 1.0;
            LambdaG_ = 1.0;
        }
        else if (VTSpeciesDependent_ and VTCollidingPartner_)
        {
            if (dict2T.subDict("strathVTCoefficients").isDict(name1+"_"+name2))
            {
                preAij = readScalar(dict2T.subDict("strathVTCoefficients").subDict(name1+"_"+name2).lookup("preAij"));
                preMij = readScalar(dict2T.subDict("strathVTCoefficients").subDict(name1+"_"+name2).lookup("preMij"));
                offset_ = readScalar(dict2T.subDict("strathVTCoefficients").subDict(name1+"_"+name2).lookup("offset"));
                sigma1_ = readScalar(dict2T.subDict("strathVTCoefficients").subDict(name1+"_"+name2).lookup("sigma1"));
                sigma2_ = readScalar(dict2T.subDict("strathVTCoefficients").subDict(name1+"_"+name2).lookup("sigma2"));
                LambdaD_ = readScalar(dict2T.subDict("strathVTCoefficients").subDict(name1+"_"+name2).lookup("LambdaD"));
                LambdaE_ = readScalar(dict2T.subDict("strathVTCoefficients").subDict(name1+"_"+name2).lookup("LambdaE"));
                LambdaG_ = readScalar(dict2T.subDict("strathVTCoefficients").subDict(name1+"_"+name2).lookup("LambdaG"));
            }
            else if (dict2T.subDict("strathVTCoefficients").isDict(name2+"_"+name1))
            {
                preAij = readScalar(dict2T.subDict("strathVTCoefficients").subDict(name2+"_"+name1).lookup("preAij"));
                preMij = readScalar(dict2T.subDict("strathVTCoefficients").subDict(name2+"_"+name1).lookup("preMij"));
                offset_ = readScalar(dict2T.subDict("strathVTCoefficients").subDict(name2+"_"+name1).lookup("offset"));
                sigma1_ = readScalar(dict2T.subDict("strathVTCoefficients").subDict(name2+"_"+name1).lookup("sigma1"));
                sigma2_ = readScalar(dict2T.subDict("strathVTCoefficients").subDict(name2+"_"+name1).lookup("sigma2"));
                LambdaD_ = readScalar(dict2T.subDict("strathVTCoefficients").subDict(name2+"_"+name1).lookup("LambdaD"));
                LambdaE_ = readScalar(dict2T.subDict("strathVTCoefficients").subDict(name2+"_"+name1).lookup("LambdaE"));
                LambdaG_ = readScalar(dict2T.subDict("strathVTCoefficients").subDict(name2+"_"+name1).lookup("LambdaG"));
            }
            else if (dict2T.subDict("strathVTCoefficients").isDict(name1))
            {
                preAij = readScalar(dict2T.subDict("strathVTCoefficients").subDict(name1).lookup("preAij"));
                preMij = readScalar(dict2T.subDict("strathVTCoefficients").subDict(name1).lookup("preMij"));
                offset_ = readScalar(dict2T.subDict("strathVTCoefficients").subDict(name1).lookup("offset"));
                sigma1_ = readScalar(dict2T.subDict("strathVTCoefficients").subDict(name1).lookup("sigma1"));
                sigma2_ = readScalar(dict2T.subDict("strathVTCoefficients").subDict(name1).lookup("sigma2"));
                LambdaD_ = readScalar(dict2T.subDict("strathVTCoefficients").subDict(name1).lookup("LambdaD"));
                LambdaE_ = readScalar(dict2T.subDict("strathVTCoefficients").subDict(name1).lookup("LambdaE"));
                LambdaG_ = readScalar(dict2T.subDict("strathVTCoefficients").subDict(name1).lookup("LambdaG"));
            }
            else
            {
                preAij = readScalar(dict2T.subDict("strathVTCoefficients").subDict("allSpecies").lookup("preAij"));
                preMij = readScalar(dict2T.subDict("strathVTCoefficients").subDict("allSpecies").lookup("preMij"));
                offset_ = readScalar(dict2T.subDict("strathVTCoefficients").subDict("allSpecies").lookup("offset"));
                sigma1_ = readScalar(dict2T.subDict("strathVTCoefficients").subDict("allSpecies").lookup("sigma1"));
                sigma2_ = readScalar(dict2T.subDict("strathVTCoefficients").subDict("allSpecies").lookup("sigma2"));
                LambdaD_ = readScalar(dict2T.subDict("strathVTCoefficients").subDict("allSpecies").lookup("LambdaD"));
                LambdaE_ = readScalar(dict2T.subDict("strathVTCoefficients").subDict("allSpecies").lookup("LambdaE"));
                LambdaG_ = readScalar(dict2T.subDict("strathVTCoefficients").subDict("allSpecies").lookup("LambdaG"));
            }
        }
        else if (VTSpeciesDependent_ and dict2T.subDict("strathVTCoefficients").isDict(name1))
        {
            preAij = readScalar(dict2T.subDict("strathVTCoefficients").subDict(name1).lookup("preAij"));
            preMij = readScalar(dict2T.subDict("strathVTCoefficients").subDict(name1).lookup("preMij"));
            offset_ = readScalar(dict2T.subDict("strathVTCoefficients").subDict(name1).lookup("offset"));
            sigma1_ = readScalar(dict2T.subDict("strathVTCoefficients").subDict(name1).lookup("sigma1"));
            sigma2_ = readScalar(dict2T.subDict("strathVTCoefficients").subDict(name1).lookup("sigma2"));
            LambdaD_ = readScalar(dict2T.subDict("strathVTCoefficients").subDict(name1).lookup("LambdaD"));
            LambdaE_ = readScalar(dict2T.subDict("strathVTCoefficients").subDict(name1).lookup("LambdaE"));
            LambdaG_ = readScalar(dict2T.subDict("strathVTCoefficients").subDict(name1).lookup("LambdaG"));
        }
        else
        {
            preAij = readScalar(dict2T.subDict("strathVTCoefficients").subDict("allSpecies").lookup("preAij"));
            preMij = readScalar(dict2T.subDict("strathVTCoefficients").subDict("allSpecies").lookup("preMij"));
            offset_ = readScalar(dict2T.subDict("strathVTCoefficients").subDict("allSpecies").lookup("offset"));
            sigma1_ = readScalar(dict2T.subDict("strathVTCoefficients").subDict("allSpecies").lookup("sigma1"));
            sigma2_ = readScalar(dict2T.subDict("strathVTCoefficients").subDict("allSpecies").lookup("sigma2"));
            LambdaD_ = readScalar(dict2T.subDict("strathVTCoefficients").subDict("allSpecies").lookup("LambdaD"));
            LambdaE_ = readScalar(dict2T.subDict("strathVTCoefficients").subDict("allSpecies").lookup("LambdaE"));
            LambdaG_ = readScalar(dict2T.subDict("strathVTCoefficients").subDict("allSpecies").lookup("LambdaG"));
        }

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
            sigma1_ = 1.0e-21;
            sigma2_ = 5.0e4;

            LambdaD_ = 0.0;
            LambdaE_ = 1.0;
            LambdaG_ = 1.0;
        }
        else if (VTSpeciesDependent_ and VTCollidingPartner_)
        {
            if (dict2T.subDict("strathVTCoefficients").isDict(name1+"_"+name2))
            {
                A12_= readScalar(dict2T.subDict("strathVTCoefficients").subDict(name1+"_"+name2).lookup("Aij"));
                B12_= readScalar(dict2T.subDict("strathVTCoefficients").subDict(name1+"_"+name2).lookup("Bij"));
                offset_ = readScalar(dict2T.subDict("strathVTCoefficients").subDict(name1+"_"+name2).lookup("offset"));
                sigma1_ = readScalar(dict2T.subDict("strathVTCoefficients").subDict(name1+"_"+name2).lookup("sigma1"));
                sigma2_ = readScalar(dict2T.subDict("strathVTCoefficients").subDict(name1+"_"+name2).lookup("sigma2"));
                LambdaD_ = readScalar(dict2T.subDict("strathVTCoefficients").subDict(name1+"_"+name2).lookup("LambdaD"));
                LambdaE_ = readScalar(dict2T.subDict("strathVTCoefficients").subDict(name1+"_"+name2).lookup("LambdaE"));
                LambdaG_ = readScalar(dict2T.subDict("strathVTCoefficients").subDict(name1+"_"+name2).lookup("LambdaG"));
            }
            else if (dict2T.subDict("strathVTCoefficients").isDict(name2+"_"+name1))
            {
                A12_= readScalar(dict2T.subDict("strathVTCoefficients").subDict(name2+"_"+name1).lookup("Aij"));
                B12_= readScalar(dict2T.subDict("strathVTCoefficients").subDict(name2+"_"+name1).lookup("Bij"));
                offset_ = readScalar(dict2T.subDict("strathVTCoefficients").subDict(name2+"_"+name1).lookup("offset"));
                sigma1_ = readScalar(dict2T.subDict("strathVTCoefficients").subDict(name2+"_"+name1).lookup("sigma1"));
                sigma2_ = readScalar(dict2T.subDict("strathVTCoefficients").subDict(name2+"_"+name1).lookup("sigma2"));
                LambdaD_ = readScalar(dict2T.subDict("strathVTCoefficients").subDict(name2+"_"+name1).lookup("LambdaD"));
                LambdaE_ = readScalar(dict2T.subDict("strathVTCoefficients").subDict(name2+"_"+name1).lookup("LambdaE"));
                LambdaG_ = readScalar(dict2T.subDict("strathVTCoefficients").subDict(name2+"_"+name1).lookup("LambdaG"));
            }
            else if (dict2T.subDict("strathVTCoefficients").isDict(name1))
            {
                A12_= readScalar(dict2T.subDict("strathVTCoefficients").subDict(name1).lookup("Aij"));
                B12_= readScalar(dict2T.subDict("strathVTCoefficients").subDict(name1).lookup("Bij"));
                offset_ = readScalar(dict2T.subDict("strathVTCoefficients").subDict(name1).lookup("offset"));
                sigma1_ = readScalar(dict2T.subDict("strathVTCoefficients").subDict(name1).lookup("sigma1"));
                sigma2_ = readScalar(dict2T.subDict("strathVTCoefficients").subDict(name1).lookup("sigma2"));
                LambdaD_ = readScalar(dict2T.subDict("strathVTCoefficients").subDict(name1).lookup("LambdaD"));
                LambdaE_ = readScalar(dict2T.subDict("strathVTCoefficients").subDict(name1).lookup("LambdaE"));
                LambdaG_ = readScalar(dict2T.subDict("strathVTCoefficients").subDict(name1).lookup("LambdaG"));
            }
            else
            {
                A12_= readScalar(dict2T.subDict("strathVTCoefficients").subDict("allSpecies").lookup("Aij"));
                B12_= readScalar(dict2T.subDict("strathVTCoefficients").subDict("allSpecies").lookup("Bij"));
                offset_ = readScalar(dict2T.subDict("strathVTCoefficients").subDict("allSpecies").lookup("offset"));
                sigma1_ = readScalar(dict2T.subDict("strathVTCoefficients").subDict("allSpecies").lookup("sigma1"));
                sigma2_ = readScalar(dict2T.subDict("strathVTCoefficients").subDict("allSpecies").lookup("sigma2"));
                LambdaD_ = readScalar(dict2T.subDict("strathVTCoefficients").subDict("allSpecies").lookup("LambdaD"));
                LambdaE_ = readScalar(dict2T.subDict("strathVTCoefficients").subDict("allSpecies").lookup("LambdaE"));
                LambdaG_ = readScalar(dict2T.subDict("strathVTCoefficients").subDict("allSpecies").lookup("LambdaG"));
            }
        }
        else if (VTSpeciesDependent_ and dict2T.subDict("strathVTCoefficients").isDict(name1))
        {
            A12_= readScalar(dict2T.subDict("strathVTCoefficients").subDict(name1).lookup("Aij"));
            B12_= readScalar(dict2T.subDict("strathVTCoefficients").subDict(name1).lookup("Bij"));
            offset_ = readScalar(dict2T.subDict("strathVTCoefficients").subDict(name1).lookup("offset"));
            sigma1_ = readScalar(dict2T.subDict("strathVTCoefficients").subDict(name1).lookup("sigma1"));
            sigma2_ = readScalar(dict2T.subDict("strathVTCoefficients").subDict(name1).lookup("sigma2"));
            LambdaD_ = readScalar(dict2T.subDict("strathVTCoefficients").subDict(name1).lookup("LambdaD"));
            LambdaE_ = readScalar(dict2T.subDict("strathVTCoefficients").subDict(name1).lookup("LambdaE"));
            LambdaG_ = readScalar(dict2T.subDict("strathVTCoefficients").subDict(name1).lookup("LambdaG"));
        }
        else
        {
            A12_= readScalar(dict2T.subDict("strathVTCoefficients").subDict("allSpecies").lookup("Aij"));
            B12_= readScalar(dict2T.subDict("strathVTCoefficients").subDict("allSpecies").lookup("Bij"));
            offset_ = readScalar(dict2T.subDict("strathVTCoefficients").subDict("allSpecies").lookup("offset"));
            sigma1_ = readScalar(dict2T.subDict("strathVTCoefficients").subDict("allSpecies").lookup("sigma1"));
            sigma2_ = readScalar(dict2T.subDict("strathVTCoefficients").subDict("allSpecies").lookup("sigma2"));
            LambdaD_ = readScalar(dict2T.subDict("strathVTCoefficients").subDict("allSpecies").lookup("LambdaD"));
            LambdaE_ = readScalar(dict2T.subDict("strathVTCoefficients").subDict("allSpecies").lookup("LambdaE"));
            LambdaG_ = readScalar(dict2T.subDict("strathVTCoefficients").subDict("allSpecies").lookup("LambdaG"));
        }
    }

}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::VTRelaxationModels::strathVT::tauVT() const
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
            1.01325e5 / this->p_[celli] * exp(A12_*(pow(this->T_[celli], -1.0/3.0) - B12_ + LambdaD_*pow(this->Tv_[species1_][celli]/THETA1_, LambdaE_)) - offset_)
          + LambdaG_/(sqrt(8.0*constant::physicoChemical::R.value()*1000.0*this->T_[celli]/
              (constant::mathematical::pi*W1_)) * sigma1_*pow(sigma2_/this->T_[celli], 2.0) *max(this->nD_[species1_][celli], Foam::SMALL));
    }


    forAll(this->T_.boundaryField(), patchi)
    {
        const fvPatchScalarField& pTt = this->T_.boundaryField()[patchi];
        const fvPatchScalarField& pp = this->p_.boundaryField()[patchi];
        const fvPatchScalarField& pTv = this->Tv_[species1_].boundaryField()[patchi];
        const fvPatchScalarField& pnD = this->nD_[species1_].boundaryField()[patchi];
        fvPatchScalarField& ptauVT = tauVT.boundaryFieldRef()[patchi];

        forAll(pTt, facei)
        {
            ptauVT[facei] =
            1.01325e5 / pp[facei] * exp(A12_*(pow(pTt[facei], -1.0/3.0) - B12_ + LambdaD_*pow(pTv[facei]/THETA1_, LambdaE_)) - offset_)
          + LambdaG_/(sqrt(8.0*constant::physicoChemical::R.value()*1000.0*pTt[facei]/
              (constant::mathematical::pi*W1_)) * sigma1_*pow(sigma2_/pTt[facei], 2.0) * max(pnD[facei], Foam::SMALL));
        }
    }

    return ttauVT;
}


Foam::tmp<Foam::scalarField> Foam::VTRelaxationModels::strathVT::tauVT
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
            1.01325e5 / p[facei] * exp(A12_*(pow(Tt[facei], -1.0/3.0) - B12_ + LambdaD_*pow(Tv[species1_][facei]/THETA1_, LambdaE_)) - offset_)
          + LambdaG_/(sqrt(8.0*constant::physicoChemical::R.value()*1000.0*Tt[facei]/
              (constant::mathematical::pi*W1_)) * sigma1_*pow(sigma2_/Tt[facei], 2.0) * max(nD[species1_][facei], Foam::SMALL));
    }

    return ttauVT;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
