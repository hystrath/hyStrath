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

#include "makeMultiSpeciesTransportModel.H"

#include "noSpeciesDiffusion.H"
#include "Fick.H"
#include "SCEBD.H"
#include "LewisNumber.H"
#include "modifiedLewisNumber.H"

#include "thermoPhysics2Types.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/* * * * * * * * * * * * * * * private static data * * * * * * * * * * * * * */

// Based on sensible enthalpy

makeMultiSpeciesTransportModel
(
    noSpeciesDiffusion,
    demConstGasHThermoPhysicsH2TGD
);
makeMultiSpeciesTransportModel(Fick, demConstGasHThermoPhysicsH2TGD);
makeMultiSpeciesTransportModel(SCEBD, demConstGasHThermoPhysicsH2TGD);
makeMultiSpeciesTransportModel(LewisNumber, demConstGasHThermoPhysicsH2TGD);
makeMultiSpeciesTransportModel
(
    modifiedLewisNumber,
    demConstGasHThermoPhysicsH2TGD
);

makeMultiSpeciesTransportModel(noSpeciesDiffusion, demGasHThermoPhysicsH2TGD);
makeMultiSpeciesTransportModel(Fick, demGasHThermoPhysicsH2TGD);
makeMultiSpeciesTransportModel(SCEBD, demGasHThermoPhysicsH2TGD);
makeMultiSpeciesTransportModel(LewisNumber, demGasHThermoPhysicsH2TGD);
makeMultiSpeciesTransportModel(modifiedLewisNumber, demGasHThermoPhysicsH2TGD);

makeMultiSpeciesTransportModel(noSpeciesDiffusion, demBEGasHThermoPhysicsH2TGD);
makeMultiSpeciesTransportModel(Fick, demBEGasHThermoPhysicsH2TGD);
makeMultiSpeciesTransportModel(SCEBD, demBEGasHThermoPhysicsH2TGD);
makeMultiSpeciesTransportModel(LewisNumber, demBEGasHThermoPhysicsH2TGD);
makeMultiSpeciesTransportModel
(
    modifiedLewisNumber,
    demBEGasHThermoPhysicsH2TGD
);

makeMultiSpeciesTransportModel
(
    noSpeciesDiffusion,
    demPLEGasHThermoPhysicsH2TGD
);
makeMultiSpeciesTransportModel(Fick, demPLEGasHThermoPhysicsH2TGD);
makeMultiSpeciesTransportModel(SCEBD, demPLEGasHThermoPhysicsH2TGD);
makeMultiSpeciesTransportModel(LewisNumber, demPLEGasHThermoPhysicsH2TGD);
makeMultiSpeciesTransportModel
(
    modifiedLewisNumber,
    demPLEGasHThermoPhysicsH2TGD
);

makeMultiSpeciesTransportModel
(
    noSpeciesDiffusion,
    demCEAGasHThermoPhysicsH2TGD
);
makeMultiSpeciesTransportModel(Fick, demCEAGasHThermoPhysicsH2TGD);
makeMultiSpeciesTransportModel(SCEBD, demCEAGasHThermoPhysicsH2TGD);
makeMultiSpeciesTransportModel(LewisNumber, demCEAGasHThermoPhysicsH2TGD);
makeMultiSpeciesTransportModel
(
    modifiedLewisNumber,
    demCEAGasHThermoPhysicsH2TGD
);

// Based on sensible internal energy

makeMultiSpeciesTransportModel
(
    noSpeciesDiffusion,
    demConstGasEThermoPhysicsH2TGD
);
makeMultiSpeciesTransportModel(Fick, demConstGasEThermoPhysicsH2TGD);
makeMultiSpeciesTransportModel(SCEBD, demConstGasEThermoPhysicsH2TGD);
makeMultiSpeciesTransportModel(LewisNumber, demConstGasEThermoPhysicsH2TGD);
makeMultiSpeciesTransportModel
(
    modifiedLewisNumber,
    demConstGasEThermoPhysicsH2TGD
);

makeMultiSpeciesTransportModel(noSpeciesDiffusion, demGasEThermoPhysicsH2TGD);
makeMultiSpeciesTransportModel(Fick, demGasEThermoPhysicsH2TGD);
makeMultiSpeciesTransportModel(SCEBD, demGasEThermoPhysicsH2TGD);
makeMultiSpeciesTransportModel(LewisNumber, demGasEThermoPhysicsH2TGD);
makeMultiSpeciesTransportModel(modifiedLewisNumber, demGasEThermoPhysicsH2TGD);

makeMultiSpeciesTransportModel(noSpeciesDiffusion, demBEGasEThermoPhysicsH2TGD);
makeMultiSpeciesTransportModel(Fick, demBEGasEThermoPhysicsH2TGD);
makeMultiSpeciesTransportModel(SCEBD, demBEGasEThermoPhysicsH2TGD);
makeMultiSpeciesTransportModel(LewisNumber, demBEGasEThermoPhysicsH2TGD);
makeMultiSpeciesTransportModel
(
    modifiedLewisNumber,
    demBEGasEThermoPhysicsH2TGD
);

makeMultiSpeciesTransportModel
(
    noSpeciesDiffusion,
    demPLEGasEThermoPhysicsH2TGD
);
makeMultiSpeciesTransportModel(Fick, demPLEGasEThermoPhysicsH2TGD);
makeMultiSpeciesTransportModel(SCEBD, demPLEGasEThermoPhysicsH2TGD);
makeMultiSpeciesTransportModel(LewisNumber, demPLEGasEThermoPhysicsH2TGD);
makeMultiSpeciesTransportModel
(
    modifiedLewisNumber,
    demPLEGasEThermoPhysicsH2TGD
);

makeMultiSpeciesTransportModel
(
    noSpeciesDiffusion,
    demCEAGasEThermoPhysicsH2TGD
);
makeMultiSpeciesTransportModel(Fick, demCEAGasEThermoPhysicsH2TGD);
makeMultiSpeciesTransportModel(SCEBD, demCEAGasEThermoPhysicsH2TGD);
makeMultiSpeciesTransportModel(LewisNumber, demCEAGasEThermoPhysicsH2TGD);
makeMultiSpeciesTransportModel
(
    modifiedLewisNumber,
    demCEAGasEThermoPhysicsH2TGD
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
