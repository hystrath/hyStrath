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

#include "SchwartzSlawskyHerzfeld.H"
#include "addToRunTimeSelectionTable.H"
#include "MillikanWhitePark.H"

#include "DynamicList.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace VTRelaxationModels
    {
        defineTypeNameAndDebug(SchwartzSlawskyHerzfeld, 0);
        addToRunTimeSelectionTable
        (
            VTRelaxationModel,
            SchwartzSlawskyHerzfeld,
            dictionary
        );
    }


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::VTRelaxationModels::SchwartzSlawskyHerzfeld::SchwartzSlawskyHerzfeld
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
    VTRelaxationModel(name1, name2, lname1, lname2, dict2T, dictThermoPhy, p, Tt, Tv, nD),

    pi(constant::mathematical::pi),
    kB(constant::physicoChemical::k.value()),
    hPlanck(constant::universal::h.value())
{
    m1_ = readScalar(dictThermoPhy.subDict(name1).subDict("specie").lookup("molWeight"))*1.0e-3/constant::physicoChemical::NA.value();
    const scalar m2 = readScalar(dictThermoPhy.subDict(name2).subDict("specie").lookup("molWeight"))*1.0e-3/constant::physicoChemical::NA.value();
    DynamicList<scalar> vibData(dictThermoPhy.subDict(name1).subDict("thermodynamics").lookup("vibrationalList"));
    TH1_ = vibData[1];
    mu12_ = reducedMass(m1_,m2);

    species1_ = lname1; species2_ = lname2;
    SHHon_ = false;

    if(dict2T.subDict("SchwartzSlawskyHerzfeldCoefficients").subDict("LennardJonesParameters").isDict(name1+"_"+name2) or dict2T.subDict("SchwartzSlawskyHerzfeldCoefficients").subDict("LennardJonesParameters").isDict(name1+"_"+name2))
    {
        SHHon_ = true;
    }

    if(SHHon_)
    {
        Tlow_ = readScalar(dict2T.subDict("SchwartzSlawskyHerzfeldCoefficients").lookup("Tlow"));
        Thigh_ = readScalar(dict2T.subDict("SchwartzSlawskyHerzfeldCoefficients").lookup("Thigh"));

        word subDictName = word::null;

        if (dict2T.subDict("SchwartzSlawskyHerzfeldCoefficients").subDict("LennardJonesParameters").isDict(name1+"_"+name2))
        {
            subDictName = name1+"_"+name2;
        }
        else
        {
            subDictName = name2+"_"+name1;
        }

        FixedList<scalar, 2> defaultList;
        defaultList[0] = 1.0e-20;
        defaultList[1] = 1.0e-20;
        sigma12_ = dict2T.subDict("SchwartzSlawskyHerzfeldCoefficients").subDict("LennardJonesParameters").subDict(subDictName).lookupOrDefault<FixedList<scalar,2> >("sigma", defaultList);
        epsilon12_ = dict2T.subDict("SchwartzSlawskyHerzfeldCoefficients").subDict("LennardJonesParameters").subDict(subDictName).lookupOrDefault<FixedList<scalar,2> >("epsilonBykB", defaultList);

        forAll(sigma12_, i)
        {
            sigma12_[i] *= 1.0e-10;
            epsilon12_[i] *= kB;
        }

        /*scalar polarizability2(0.0);
        if (dict2T.subDict("SchwartzSlawskyHerzfeldCoefficients").subDict("speciesData").isDict(name2))
        {
            polarizability2 = readScalar(dict2T.subDict("SchwartzSlawskyHerzfeldCoefficients").subDict("speciesData").subDict(name2).lookup("polarizability"));
        }
        else
        {
            FatalErrorIn("Foam::VTRelaxationModels::SchwartzSlawskyHerzfeld::SchwartzSlawskyHerzfeld")
                        << "Some coefficients are missing." << nl;
            FatalError<< exit(FatalError);
        }*/ // DEACTIVATED VINCENT 10/03/2016

        if (dict2T.subDict("SchwartzSlawskyHerzfeldCoefficients").subDict("speciesData").isDict(name1))
        {
            rp1_ = 1.0e-10*readScalar(dict2T.subDict("SchwartzSlawskyHerzfeldCoefficients").subDict("speciesData").subDict(name1).lookup("LVLInteratomicDistance"));
        }
        else
        {
            FatalErrorIn("Foam::VTRelaxationModels::SchwartzSlawskyHerzfeld::SchwartzSlawskyHerzfeld")
                        << "Some SSH coefficients are missing." << nl;
            FatalError<< exit(FatalError);
        }

        /*forAll(sigma12_, i)
        {
            sigma12_[i] = sigmaLJmixture(name1, m1_, sigma1[i], epsilon1[i], name2, sigma2[i], epsilon2[i], polarizability2);
            epsilon12_[i] = epsilonLJmixture(m1_, sigma1[i], epsilon1[i], sigma2[i], epsilon2[i], polarizability2);
        }*/ // DEACTIVATED VINCENT 10/03/2016

        matchingSpeciesIndices_ = similarSpecies(lname1, lname2);
        fms_ = fms(name1);
    }
    else
    {
        W1_ = readScalar(dictThermoPhy.subDict(name1).subDict("specie").lookup("molWeight"));

        word subDictName = word::null;

        if (not VTFullCoeffsForm_)
        {
            const scalar W2 = readScalar(dictThermoPhy.subDict(name2).subDict("specie").lookup("molWeight"));
            DynamicList<scalar> vibData(dictThermoPhy.subDict(name1).subDict("thermodynamics").lookup("vibrationalList"));
            const scalar TH1 = vibData[1];

            scalar W12 = (W1_ * W2) / (W1_ + W2);
            A12_ = sqrt(W12) * pow(TH1, 4.0/3.0);
            B12_ = pow025(W12);
            scalar preAij = 0.0;
            scalar preMij = 0.0;

            if (not VTOverwriteDefault_)
            {
                preAij  = 1.16e-3;
                preMij  = 0.015;
                offset_ = 18.42;
                sigma1_ = 1.0e-21;
                sigma2_ = 5.0e4;
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
            }
        }
    }
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::VTRelaxationModels::SchwartzSlawskyHerzfeld::tauVT() const
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

    if(SHHon_)
    {
        forAll(this->T_, celli)
        {
            tauVT[celli] = 1.0 /
                (
                    (1.0-exp(-TH1_/this->T_[celli]))
                  * P10sr(this->T_[celli])
                  * Zcollsr(this->T_[celli], this->nD_[species1_][celli], this->nD_[species2_][celli])
                );
        }


        forAll(this->T_.boundaryField(), patchi)
        {
            const fvPatchScalarField& pTt = this->T_.boundaryField()[patchi];
            const fvPatchScalarField& pnD1 = this->nD_[species1_].boundaryField()[patchi];
            const fvPatchScalarField& pnD2 = this->nD_[species2_].boundaryField()[patchi];
            fvPatchScalarField& ptauVT = tauVT.boundaryFieldRef()[patchi];

            forAll(pTt, facei)
            {
                ptauVT[facei] = 1.0 /
                (
                    (1.0-exp(-TH1_/pTt[facei]))
                  * P10sr(pTt[facei])
                  * Zcollsr(pTt[facei], pnD1[facei], pnD2[facei])
                );
            }
        }
    }
    else
    {
        volScalarField nDcol = this->nD_[species1_];
        if(species1_ != species2_)
        {
            nDcol += this->nD_[species2_];
        }

        forAll(this->T_, celli)
        {
            tauVT[celli] =
                1.01325e5 / this->p_[celli] * exp(A12_*(pow(this->T_[celli], -1.0/3.0) - B12_) - offset_)
              + 1.0/(sqrt(8.0*constant::physicoChemical::R.value()*1000.0*this->T_[celli]/
                  (pi*W1_)) * sigma1_*sqr(sigma2_/this->T_[celli]) *max(nDcol[celli], Foam::SMALL));
        }


        forAll(this->T_.boundaryField(), patchi)
        {
            const fvPatchScalarField& pTt = this->T_.boundaryField()[patchi];
            const fvPatchScalarField& pp = this->p_.boundaryField()[patchi];
            const fvPatchScalarField& pnDcol = nDcol.boundaryField()[patchi];
            fvPatchScalarField& ptauVT = tauVT.boundaryFieldRef()[patchi];

            forAll(pTt, facei)
            {
                ptauVT[facei] =
                1.01325e5 / pp[facei] * exp(A12_*(pow(pTt[facei], -1.0/3.0) - B12_) - offset_)
              + 1.0/(sqrt(8.0*constant::physicoChemical::R.value()*1000.0*pTt[facei]/
                  (pi*W1_)) * sigma1_*sqr(sigma2_/pTt[facei]) * max(pnDcol[facei], Foam::SMALL));
            }
        }
    }

    return ttauVT;
}


Foam::tmp<Foam::scalarField> Foam::VTRelaxationModels::SchwartzSlawskyHerzfeld::tauVT
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

    if(SHHon_)
    {
        forAll(Tt, facei)
        {
            tauVT[facei] = 1.0 /
                (
                    (1.0-exp(-TH1_/Tt[facei]))
                  * P10sr(Tt[facei])
                  * Zcollsr(Tt[facei], nD[species1_][facei], nD[species2_][facei])
                );
        }
    }
    else
    {
        scalarField nDcol = nD[species1_];
        if(species1_ != species2_)
        {
            nDcol += nD[species2_];
        }

        forAll(Tt, facei)
        {
            tauVT[facei] =
                1.01325e5 / p[facei] * exp(A12_*(pow(Tt[facei], -1.0/3.0) - B12_) - offset_)
              + 1.0/(sqrt(8.0*constant::physicoChemical::R.value()*1000.0*Tt[facei]/
                  (pi*W1_)) * sigma1_*sqr(sigma2_/Tt[facei]) * max(nDcol[facei],Foam::SMALL));
        }
    }

    return ttauVT;
}


} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
