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

#include "makeRelaxationTimeModelVV.H"

#include "noVVEnergyTransfer.H"
#include "KnabVV.H"

#include "thermoPhysics2Types.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/* * * * * * * * * * * * * * * private static data * * * * * * * * * * * * * */

// Based on sensible enthalpy

makeRelaxationTimeModelVV(noVVEnergyTransfer, demConstGasHThermoPhysicsH2TGD);
makeRelaxationTimeModelVV(KnabVV, demConstGasHThermoPhysicsH2TGD);

makeRelaxationTimeModelVV(noVVEnergyTransfer, demGasHThermoPhysicsH2TGD);
makeRelaxationTimeModelVV(KnabVV, demGasHThermoPhysicsH2TGD);

makeRelaxationTimeModelVV(noVVEnergyTransfer, demBEGasHThermoPhysicsH2TGD);
makeRelaxationTimeModelVV(KnabVV, demBEGasHThermoPhysicsH2TGD);

makeRelaxationTimeModelVV(noVVEnergyTransfer, demPLEGasHThermoPhysicsH2TGD);
makeRelaxationTimeModelVV(KnabVV, demPLEGasHThermoPhysicsH2TGD);

makeRelaxationTimeModelVV(noVVEnergyTransfer, demCEAGasHThermoPhysicsH2TGD);
makeRelaxationTimeModelVV(KnabVV, demCEAGasHThermoPhysicsH2TGD);

// Based on sensible internal energy

makeRelaxationTimeModelVV(noVVEnergyTransfer, demConstGasEThermoPhysicsH2TGD);
makeRelaxationTimeModelVV(KnabVV, demConstGasEThermoPhysicsH2TGD);

makeRelaxationTimeModelVV(noVVEnergyTransfer, demGasEThermoPhysicsH2TGD);
makeRelaxationTimeModelVV(KnabVV, demGasEThermoPhysicsH2TGD);

makeRelaxationTimeModelVV(noVVEnergyTransfer, demBEGasEThermoPhysicsH2TGD);
makeRelaxationTimeModelVV(KnabVV, demBEGasEThermoPhysicsH2TGD);

makeRelaxationTimeModelVV(noVVEnergyTransfer, demPLEGasEThermoPhysicsH2TGD);
makeRelaxationTimeModelVV(KnabVV, demPLEGasEThermoPhysicsH2TGD);

makeRelaxationTimeModelVV(noVVEnergyTransfer, demCEAGasEThermoPhysicsH2TGD);
makeRelaxationTimeModelVV(KnabVV, demCEAGasEThermoPhysicsH2TGD);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
