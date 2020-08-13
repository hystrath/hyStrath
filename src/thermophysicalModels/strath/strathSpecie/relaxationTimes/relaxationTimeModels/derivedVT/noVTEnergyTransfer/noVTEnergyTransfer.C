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

#include "noVTEnergyTransfer.H"
#include "fvm.H"

// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

template<class ThermoType>
void Foam::noVTEnergyTransfer<ThermoType>::updateCoefficients()
{
    forAll(species(), i)
    {
        tauVT_[i] = dimensionedScalar("GREAT", dimTime, Foam::GREAT);
    }

    /*forAll(species(), i) // TODO ABORTIVE WORK
    {
        forAll(tauVTmode_[i], m)
        {
            tauVTmode_[i][m] = dimensionedScalar("GREAT", dimTime, Foam::GREAT);
        }
    }*/
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ThermoType>
Foam::noVTEnergyTransfer<ThermoType>::noVTEnergyTransfer
(
    rho2ReactionThermo& thermo,
    const compressibleTurbulenceModel& turbulence
)
:
    relaxationTimeModel(thermo, turbulence),

    speciesThermo_
    (
        dynamic_cast<const multi2ComponentMixture<ThermoType>&>
            (this->thermo_).speciesData()
    )
{
    tauVT_.setSize(solvedVibEqSpecies().size());
    //tauVTmode_.setSize(species().size()); // TODO ABORTIVE WORK

    forAll(tauVT_, speciei)
    {
        tauVT_.set
        (
            speciei,
            new volScalarField
            (
                IOobject
                (
                    "tauVT_" + species()[speciei],
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh_,
                dimensionedScalar("tauVT", dimTime, 0.0)
            )
        );
    }

    /*forAll(tauVTmode_, speciei) // TODO ABORTIVE WORK
    {
        tauVTmode_.set
        (
            speciei,
            new PtrList<volScalarField>(thermo_.composition().noVibrationalTemp(speciei))
        );
    }

    forAll(tauVTmode_, speciei)
    {
      forAll(tauVTmode_[speciei], vibMode)
      {
        tauVTmode_[speciei].set
        (
            vibMode,
            new volScalarField
            (
                IOobject
                (
                    "tauVT_" + species()[speciei] + "." + word(vibMode+1),
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh_,
                dimensionedScalar("tauVT", dimTime, 0.0)
            )
        );
      }
    }*/
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class ThermoType>
void Foam::noVTEnergyTransfer<ThermoType>::correct()
{}


template<class ThermoType>
bool Foam::noVTEnergyTransfer<ThermoType>::read()
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
