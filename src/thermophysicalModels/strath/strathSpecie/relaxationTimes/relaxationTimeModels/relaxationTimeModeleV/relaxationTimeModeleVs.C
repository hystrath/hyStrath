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

#include "makeRelaxationTimeModeleV.H"

#include "noeVEnergyTransfer.H"
#include "LeeLandauTellereV.H"

#include "thermoPhysics2Types.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/* * * * * * * * * * * * * * * private static data * * * * * * * * * * * * * */

// Based on sensible enthalpy

makeRelaxationTimeModeleV(noeVEnergyTransfer, demConstGasHThermoPhysicsH2TGD);
makeRelaxationTimeModeleV(LeeLandauTellereV, demConstGasHThermoPhysicsH2TGD);

makeRelaxationTimeModeleV(noeVEnergyTransfer, demGasHThermoPhysicsH2TGD);
makeRelaxationTimeModeleV(LeeLandauTellereV, demGasHThermoPhysicsH2TGD);

makeRelaxationTimeModeleV(noeVEnergyTransfer, demBEGasHThermoPhysicsH2TGD);
makeRelaxationTimeModeleV(LeeLandauTellereV, demBEGasHThermoPhysicsH2TGD);

makeRelaxationTimeModeleV(noeVEnergyTransfer, demPLEGasHThermoPhysicsH2TGD);
makeRelaxationTimeModeleV(LeeLandauTellereV, demPLEGasHThermoPhysicsH2TGD);

makeRelaxationTimeModeleV(noeVEnergyTransfer, demCEAGasHThermoPhysicsH2TGD);
makeRelaxationTimeModeleV(LeeLandauTellereV, demCEAGasHThermoPhysicsH2TGD);

// Based on sensible internal energy

makeRelaxationTimeModeleV(noeVEnergyTransfer, demConstGasEThermoPhysicsH2TGD);
makeRelaxationTimeModeleV(LeeLandauTellereV, demConstGasEThermoPhysicsH2TGD);

makeRelaxationTimeModeleV(noeVEnergyTransfer, demGasEThermoPhysicsH2TGD);
makeRelaxationTimeModeleV(LeeLandauTellereV, demGasEThermoPhysicsH2TGD);

makeRelaxationTimeModeleV(noeVEnergyTransfer, demBEGasEThermoPhysicsH2TGD);
makeRelaxationTimeModeleV(LeeLandauTellereV, demBEGasEThermoPhysicsH2TGD);

makeRelaxationTimeModeleV(noeVEnergyTransfer, demPLEGasEThermoPhysicsH2TGD);
makeRelaxationTimeModeleV(LeeLandauTellereV, demPLEGasEThermoPhysicsH2TGD);

makeRelaxationTimeModeleV(noeVEnergyTransfer, demCEAGasEThermoPhysicsH2TGD);
makeRelaxationTimeModeleV(LeeLandauTellereV, demCEAGasEThermoPhysicsH2TGD);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
