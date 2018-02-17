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

#include "makeRelaxationTimeModel.H"

#include "noVTEnergyTransfer.H"
#include "LandauTellerVT.H"
//#include "LarsenBorgnakkeVT.H"

#include "thermoPhysics2Types.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/* * * * * * * * * * * * * * * private static data * * * * * * * * * * * * * */

// Based on sensible enthalpy

makeRelaxationTimeModel(noVTEnergyTransfer, demConstGasHThermoPhysicsH2TGD);
makeRelaxationTimeModel(LandauTellerVT, demConstGasHThermoPhysicsH2TGD);
//makeRelaxationTimeModel(LarsenBorgnakkeVT, demConstGasHThermoPhysicsH2TGD);

makeRelaxationTimeModel(noVTEnergyTransfer, demGasHThermoPhysicsH2TGD);
makeRelaxationTimeModel(LandauTellerVT, demGasHThermoPhysicsH2TGD);
//makeRelaxationTimeModel(LarsenBorgnakkeVT, demGasHThermoPhysicsH2TGD);

makeRelaxationTimeModel(noVTEnergyTransfer, demBEGasHThermoPhysicsH2TGD);
makeRelaxationTimeModel(LandauTellerVT, demBEGasHThermoPhysicsH2TGD);
//makeRelaxationTimeModel(LarsenBorgnakkeVT, demBEGasHThermoPhysicsH2TGD);

makeRelaxationTimeModel(noVTEnergyTransfer, demPLEGasHThermoPhysicsH2TGD);
makeRelaxationTimeModel(LandauTellerVT, demPLEGasHThermoPhysicsH2TGD);
//makeRelaxationTimeModel(LarsenBorgnakkeVT, demPLEGasHThermoPhysicsH2TGD);

makeRelaxationTimeModel(noVTEnergyTransfer, demCEAGasHThermoPhysicsH2TGD);
makeRelaxationTimeModel(LandauTellerVT, demCEAGasHThermoPhysicsH2TGD);
//makeRelaxationTimeModel(LarsenBorgnakkeVT, demCEAGasHThermoPhysicsH2TGD);

// Based on sensible internal energy

makeRelaxationTimeModel(noVTEnergyTransfer, demConstGasEThermoPhysicsH2TGD);
makeRelaxationTimeModel(LandauTellerVT, demConstGasEThermoPhysicsH2TGD);
//makeRelaxationTimeModel(LarsenBorgnakkeVT, demConstGasEThermoPhysicsH2TGD);

makeRelaxationTimeModel(noVTEnergyTransfer, demGasEThermoPhysicsH2TGD);
makeRelaxationTimeModel(LandauTellerVT, demGasEThermoPhysicsH2TGD);
//makeRelaxationTimeModel(LarsenBorgnakkeVT, demGasEThermoPhysicsH2TGD);

makeRelaxationTimeModel(noVTEnergyTransfer, demBEGasEThermoPhysicsH2TGD);
makeRelaxationTimeModel(LandauTellerVT, demBEGasEThermoPhysicsH2TGD);
//makeRelaxationTimeModel(LarsenBorgnakkeVT, demBEGasEThermoPhysicsH2TGD);

makeRelaxationTimeModel(noVTEnergyTransfer, demPLEGasEThermoPhysicsH2TGD);
makeRelaxationTimeModel(LandauTellerVT, demPLEGasEThermoPhysicsH2TGD);
//makeRelaxationTimeModel(LarsenBorgnakkeVT, demPLEGasEThermoPhysicsH2TGD);

makeRelaxationTimeModel(noVTEnergyTransfer, demCEAGasEThermoPhysicsH2TGD);
makeRelaxationTimeModel(LandauTellerVT, demCEAGasEThermoPhysicsH2TGD);
//makeRelaxationTimeModel(LarsenBorgnakkeVT, demCEAGasEThermoPhysicsH2TGD);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
