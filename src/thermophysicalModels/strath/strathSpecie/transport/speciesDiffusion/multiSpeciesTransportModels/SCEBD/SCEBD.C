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

#include "SCEBD.H"
//#include <ctime>

// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

template<class ThermoType>
void Foam::SCEBD<ThermoType>::updateCoefficients()
{
    DijModel_().update();

    forAll(species(), speciei)
    {
        volScalarField tmpSum = 0 / Dij(0,0);

        const volScalarField omegai = thermo_.composition().pD(speciei)/sqrt(W(speciei));
        volScalarField omega = omegai;

        forAll(species(), speciej)
        {
            if (speciej != speciei and thermo_.composition().particleType(speciej) != 0)
            {
                tmpSum += thermo_.composition().X(speciej) / Dij(speciei, speciej);

                omega += thermo_.composition().pD(speciej)/sqrt(W(speciej));
            }
        }

        D_[speciei] = thermo_.rho()*(1.0 - omegai/omega)
            / (tmpSum + dimensionedScalar("VSMALL", dimTime/dimArea, Foam::VSMALL));

        const volScalarField& Xi = thermo_.composition().X(speciei);

        forAll(D_[speciei], celli)
        {
            if (1.0 - Xi[celli] < miniXs_)
            {
                D_[speciei][celli] = 0;
            }
        }

        forAll(D_[speciei].boundaryField(), patchi)
        {
            forAll(D_[speciei].boundaryField()[patchi], facei)
            {
                if (1.0 - Xi.boundaryField()[patchi][facei] < miniXs_)
                {
                    D_[speciei].boundaryFieldRef()[patchi][facei] = 0;
                }
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ThermoType>
Foam::SCEBD<ThermoType>::SCEBD
(
    rho2ReactionThermo& thermo,
    const compressible::turbulenceModel& turbulence
)
:
    multiSpeciesTransportModel(thermo, turbulence),

    speciesThermo_
    (
        dynamic_cast<const multi2ComponentMixture<ThermoType>&>
            (this->thermo_).speciesData()
    ),

    miniXs_(1.0e-12)
{
    D_.setSize(species().size());

    forAll(species(), speciei)
    {
        D_.set
        (
            speciei,
            new volScalarField
            (
                IOobject
                (
                    "rhoD_" + species()[speciei],
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh_,
                dimensionedScalar("D", dimMass/dimLength/dimTime, 0.0)
            )
        );
    }
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class ThermoType>
void Foam::SCEBD<ThermoType>::correct()
{
    updateCoefficients();

    if(addPressureGradientTerm_)
    {
        pressureGradientContributionToSpeciesMassFlux();
    }

    if(addTemperatureGradientTerm_)
    {
        temperatureGradientContributionToSpeciesMassFlux();
    }

    forAll(species(), speciei)
    {
        calculateJ(speciei);
    }

    calculateSumDiffusionFluxes();
}


template<class ThermoType>
bool Foam::SCEBD<ThermoType>::read()
{
    if (regIOobject::read())
    {
        return true;
    }
    else
    {
        return false;
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
