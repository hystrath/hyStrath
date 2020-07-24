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

#include "modifiedLewisNumber.H"

// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

template<class ThermoType>
void Foam::modifiedLewisNumber<ThermoType>::updateCoefficients()
{
    forAll(this->D_, speciei)
    {
        volScalarField constantrhoD
        (
            IOobject
            (
                "constantrhoD",
                this->mesh_.time().timeName(),
                this->mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            this->mesh_,
            // In this case, rho*Ds = cste independent from the temperature as opposed to the LewisNumber class
            dimensionedScalar("constantrhoD", dimMass/dimLength/dimTime, Le_)
        );

        constantrhoD.primitiveFieldRef() *= this->thermo_.rho().internalField();
        constantrhoD.boundaryFieldRef() *= this->thermo_.rho().boundaryField();

        this->D_[speciei] = constantrhoD;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ThermoType>
Foam::modifiedLewisNumber<ThermoType>::modifiedLewisNumber
(
    rho2ReactionThermo& thermo,
    const compressible::turbulenceModel& turbulence
)
:
    Fick<ThermoType>(thermo, turbulence),

    Le_(readScalar(IOdictionary::subDict("transportModels")
        .subDict("diffusionModelParameters").lookup("modifiedLewisNumber")))
{}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
