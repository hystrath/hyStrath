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

#include "makeRelaxationTimeModelHE.H"

#include "noHEEnergyTransfer.H"
#include "AppletonBray.H"

#include "thermoPhysics2Types.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/* * * * * * * * * * * * * * * private static data * * * * * * * * * * * * * */

// Based on sensible enthalpy

makeRelaxationTimeModelHE(noHEEnergyTransfer, demConstGasHThermoPhysicsH2TGD);
makeRelaxationTimeModelHE(AppletonBray, demConstGasHThermoPhysicsH2TGD);

makeRelaxationTimeModelHE(noHEEnergyTransfer, demGasHThermoPhysicsH2TGD);
makeRelaxationTimeModelHE(AppletonBray, demGasHThermoPhysicsH2TGD);

makeRelaxationTimeModelHE(noHEEnergyTransfer, demBEGasHThermoPhysicsH2TGD);
makeRelaxationTimeModelHE(AppletonBray, demBEGasHThermoPhysicsH2TGD);

makeRelaxationTimeModelHE(noHEEnergyTransfer, demPLEGasHThermoPhysicsH2TGD);
makeRelaxationTimeModelHE(AppletonBray, demPLEGasHThermoPhysicsH2TGD);

makeRelaxationTimeModelHE(noHEEnergyTransfer, demCEAGasHThermoPhysicsH2TGD);
makeRelaxationTimeModelHE(AppletonBray, demCEAGasHThermoPhysicsH2TGD);

// Based on sensible internal energy

makeRelaxationTimeModelHE(noHEEnergyTransfer, demConstGasEThermoPhysicsH2TGD);
makeRelaxationTimeModelHE(AppletonBray, demConstGasEThermoPhysicsH2TGD);

makeRelaxationTimeModelHE(noHEEnergyTransfer, demGasEThermoPhysicsH2TGD);
makeRelaxationTimeModelHE(AppletonBray, demGasEThermoPhysicsH2TGD);

makeRelaxationTimeModelHE(noHEEnergyTransfer, demBEGasEThermoPhysicsH2TGD);
makeRelaxationTimeModelHE(AppletonBray, demBEGasEThermoPhysicsH2TGD);

makeRelaxationTimeModelHE(noHEEnergyTransfer, demPLEGasEThermoPhysicsH2TGD);
makeRelaxationTimeModelHE(AppletonBray, demPLEGasEThermoPhysicsH2TGD);

makeRelaxationTimeModelHE(noHEEnergyTransfer, demCEAGasEThermoPhysicsH2TGD);
makeRelaxationTimeModelHE(AppletonBray, demCEAGasEThermoPhysicsH2TGD);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
