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

#include "strathMWP.H"
#include "addToRunTimeSelectionTable.H"

#include "DynamicList.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace VTRelaxationModels
    {
        defineTypeNameAndDebug(strathMWP, 0);
        addToRunTimeSelectionTable
        (
            VTRelaxationModel,
            strathMWP,
            dictionary
        );
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::VTRelaxationModels::strathMWP::strathMWP
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
    species1_ = lname1; species2_ = lname2;
    W1_ = readScalar(dictThermoPhy.subDict(name1).subDict("specie").lookup("molWeight"));

    word subDictName = word::null;

    if (not VTFullCoeffsForm_)
    {
        const scalar W2 = readScalar(dictThermoPhy.subDict(name2).subDict("specie").lookup("molWeight"));
        DynamicList<scalar> vibData(dictThermoPhy.subDict(name1).subDict("thermodynamics").lookup("vibrationalList"));
        const scalar TH1 = vibData[1];

        scalar W12 = (W1_ * W2) / (W1_ + W2);
        A12_ = sqrt(W12) * pow(TH1, 4.0/3.0);
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
            sigma1v_ = sigma1_;
            sigma2v_ = sigma2_;
        }
        else
        {
            if (VTSpeciesDependent_ and VTCollidingPartner_)
            {
                if (dict2T.subDict("ParkCoefficients").isDict(name1+"_"+name2))
                {
                    subDictName = name1+"_"+name2;
                }
                else if (dict2T.subDict("ParkCoefficients").isDict(name2+"_"+name1))
                {
                    subDictName = name2+"_"+name1;
                }
                else if (dict2T.subDict("ParkCoefficients").isDict(name1))
                {
                    subDictName = name1;
                }
                else
                {
                    subDictName = "allSpecies";
                }
            }
            else if (VTSpeciesDependent_ and dict2T.subDict("ParkCoefficients").isDict(name1))
            {
                subDictName = name1;
            }
            else
            {
                subDictName = "allSpecies";
            }

            preAij = readScalar(dict2T.subDict("ParkCoefficients").subDict(subDictName).lookup("preAij"));
            preMij = readScalar(dict2T.subDict("ParkCoefficients").subDict(subDictName).lookup("preMij"));

            A12_ *= preAij;
            B12_ *= preMij;

            offset_ = readScalar(dict2T.subDict("ParkCoefficients").subDict(subDictName).lookup("offset"));
            sigma1_ = readScalar(dict2T.subDict("ParkCoefficients").subDict(subDictName).lookup("sigma1"));
            sigma2_ = readScalar(dict2T.subDict("ParkCoefficients").subDict(subDictName).lookup("sigma2"));
            sigma1v_ = readScalar(dict2T.subDict("ParkCoefficients").subDict(subDictName).lookup("sigma1v"));
            sigma2v_ = readScalar(dict2T.subDict("ParkCoefficients").subDict(subDictName).lookup("sigma2v"));
        }
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
            sigma1v_ = sigma1_;
            sigma2v_ = sigma2_;
        }
        else
        {
            if (VTSpeciesDependent_ and VTCollidingPartner_)
            {
                if (dict2T.subDict("ParkCoefficients").isDict(name1+"_"+name2))
                {
                    subDictName = name1+"_"+name2;
                }
                else if (dict2T.subDict("ParkCoefficients").isDict(name2+"_"+name1))
                {
                    subDictName = name2+"_"+name1;
                }
                else if (dict2T.subDict("ParkCoefficients").isDict(name1))
                {
                    subDictName = name1;
                }
                else
                {
                    subDictName = "allSpecies";
                }
            }
            else if (VTSpeciesDependent_ and dict2T.subDict("ParkCoefficients").isDict(name1))
            {
                subDictName = name1;
            }
            else
            {
                subDictName = "allSpecies";
            }

            A12_ = readScalar(dict2T.subDict("ParkCoefficients").subDict(subDictName).lookup("Aij"));
            B12_ = readScalar(dict2T.subDict("ParkCoefficients").subDict(subDictName).lookup("Bij"));

            offset_ = readScalar(dict2T.subDict("ParkCoefficients").subDict(subDictName).lookup("offset"));
            sigma1_ = readScalar(dict2T.subDict("ParkCoefficients").subDict(subDictName).lookup("sigma1"));
            sigma2_ = readScalar(dict2T.subDict("ParkCoefficients").subDict(subDictName).lookup("sigma2"));
            sigma1v_ = readScalar(dict2T.subDict("ParkCoefficients").subDict(subDictName).lookup("sigma1v"));
            sigma2v_ = readScalar(dict2T.subDict("ParkCoefficients").subDict(subDictName).lookup("sigma2v"));
        }
    }
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::VTRelaxationModels::strathMWP::tauVT() const
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
            dimensionSet(0, 0, 1, 0, 0)
        )
    );

    volScalarField& tauVT = ttauVT.ref();

    volScalarField nDcol = this->nD_[species1_];
    if(species1_ != species2_)
    {
        nDcol += this->nD_[species2_];
    }

    volScalarField myT = this->Tv_[species1_]*0;

    forAll(this->T_, celli)
    {
        myT[celli] = max(this->Tv_[species1_][celli], this->T_[celli]);
        scalar Park = 0.0;
        if(this->T_[celli] >= this->Tv_[species1_][celli])
        {
            Park = 1.0/(sqrt(8.0*constant::physicoChemical::R.value()*1000.0*this->T_[celli]/
              (constant::mathematical::pi*W1_)) * sigma1_*pow(sigma2_/this->T_[celli], 2.0) *max(nDcol[celli], Foam::SMALL));
        }
        else
        {
            Park = 1.0/(sqrt(8.0*constant::physicoChemical::R.value()*1000.0*myT[celli]/
              (constant::mathematical::pi*W1_)) * (sigma1v_-sigma1v_)*pow((sigma2_-sigma2v_)/myT[celli], 2.0) *max(nDcol[celli], Foam::SMALL));
        }

        tauVT[celli] =
            1.01325e5 / this->p_[celli] * exp(A12_*(pow(myT[celli], -1.0/3.0) - B12_) - offset_)
          + Park;
    }


    forAll(this->T_.boundaryField(), patchi)
    {
        const fvPatchScalarField& pTt = this->T_.boundaryField()[patchi];
        const fvPatchScalarField& pmyTv = this->Tv_[species1_].boundaryField()[patchi];
        const fvPatchScalarField& pp = this->p_.boundaryField()[patchi];
        const fvPatchScalarField& pnDcol = nDcol.boundaryField()[patchi];
        fvPatchScalarField& ptauVT = tauVT.boundaryFieldRef()[patchi];

        forAll(pTt, facei)
        {
            ptauVT[facei] =
            1.01325e5 / pp[facei] * exp(A12_*(pow(max(pTt[facei],pmyTv[facei]), -1.0/3.0) - B12_) - offset_)
          + 1.0/(sqrt(8.0*constant::physicoChemical::R.value()*1000.0*max(pTt[facei],pmyTv[facei])/
              (constant::mathematical::pi*W1_)) * sigma1_*pow(sigma2_/max(pTt[facei],pmyTv[facei]), 2.0) * max(pnDcol[facei], Foam::SMALL));
        }
    }

    return ttauVT;
}


Foam::tmp<Foam::scalarField> Foam::VTRelaxationModels::strathMWP::tauVT
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

    scalarField nDcol = nD[species1_];
    if(species1_ != species2_)
    {
        nDcol += nD[species2_];
    }

    forAll(Tt, facei)
    {
        tauVT[facei] =
            1.01325e5 / p[facei] * exp(A12_*(pow(max(Tt[facei], Tv[species1_][facei]), -1.0/3.0) - B12_) - offset_)
          + 1.0/(sqrt(8.0*constant::physicoChemical::R.value()*1000.0*max(Tt[facei], Tv[species1_][facei])/
              (constant::mathematical::pi*W1_)) * sigma1_*pow(sigma2_/max(Tt[facei], Tv[species1_][facei]), 2.0) * max(nDcol[facei],Foam::SMALL));
    }

    return ttauVT;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
