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

InClass
    Foam::rho2ChemistryModel

Description
    Creates chemistry model instances templated on the type of thermodynamics

\*---------------------------------------------------------------------------*/

#include "make2ChemistryModel.H"

#include "rho2ChemistryModel.H"
#include "chemistry2Model.H"

#include "thermoPhysics2Types.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    // Chemistry models based on sensibleEnthalpy

    make2ChemistryModel
    (
        chemistry2Model,
        rho2ChemistryModel,
        demConstGasHThermoPhysicsH2TGD
    );
    make2ChemistryModel
    (
        chemistry2Model,
        rho2ChemistryModel,
        demGasHThermoPhysicsH2TGD
    );

    make2ChemistryModel
    (
        chemistry2Model,
        rho2ChemistryModel,
        demBEGasHThermoPhysicsH2TGD
    );

    make2ChemistryModel
    (
        chemistry2Model,
        rho2ChemistryModel,
        demPLEGasHThermoPhysicsH2TGD
    );

    make2ChemistryModel
    (
        chemistry2Model,
        rho2ChemistryModel,
        demCEAGasHThermoPhysicsH2TGD
    );

    // Chemistry moldels based on sensibleInternalEnergy
    make2ChemistryModel
    (
        chemistry2Model,
        rho2ChemistryModel,
        demConstGasEThermoPhysicsH2TGD
    );

    make2ChemistryModel
    (
        chemistry2Model,
        rho2ChemistryModel,
        demGasEThermoPhysicsH2TGD
    );

    make2ChemistryModel
    (
        chemistry2Model,
        rho2ChemistryModel,
        demBEGasEThermoPhysicsH2TGD
    );

    make2ChemistryModel
    (
        chemistry2Model,
        rho2ChemistryModel,
        demPLEGasEThermoPhysicsH2TGD
    );

    make2ChemistryModel
    (
        chemistry2Model,
        rho2ChemistryModel,
        demCEAGasEThermoPhysicsH2TGD
    );

}

// ************************************************************************* //
