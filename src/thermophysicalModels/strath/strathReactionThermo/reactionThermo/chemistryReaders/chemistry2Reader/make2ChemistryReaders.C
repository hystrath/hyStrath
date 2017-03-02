/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

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

#include "make2ReactionThermo.H"
#include "thermoPhysics2Types.H"
#include "solidThermoPhysicsTypes.H"

#include "chemistry2Reader.H"
#include "foam2ChemistryReader.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Solid chemistry readers based on sensibleEnthalpy

make2ChemistryReader(demConstGasHThermoPhysicsH2TGD);
make2ChemistryReaderType(foam2ChemistryReader, demConstGasHThermoPhysicsH2TGD);

make2ChemistryReader(demGasHThermoPhysicsH2TGD);
make2ChemistryReaderType(foam2ChemistryReader, demGasHThermoPhysicsH2TGD);

make2ChemistryReader(demBEGasHThermoPhysicsH2TGD);
make2ChemistryReaderType(foam2ChemistryReader, demBEGasHThermoPhysicsH2TGD);

make2ChemistryReader(demPLEGasHThermoPhysicsH2TGD);
make2ChemistryReaderType(foam2ChemistryReader, demPLEGasHThermoPhysicsH2TGD);

make2ChemistryReader(demCEAGasHThermoPhysicsH2TGD);
make2ChemistryReaderType(foam2ChemistryReader, demCEAGasHThermoPhysicsH2TGD);


// Solid chemistry readers for solids based on sensibleInternalEnergy

make2ChemistryReader(demConstGasEThermoPhysicsH2TGD);
make2ChemistryReaderType(foam2ChemistryReader, demConstGasEThermoPhysicsH2TGD);

make2ChemistryReader(demGasEThermoPhysicsH2TGD);
make2ChemistryReaderType(foam2ChemistryReader, demGasEThermoPhysicsH2TGD);

make2ChemistryReader(demBEGasEThermoPhysicsH2TGD);
make2ChemistryReaderType(foam2ChemistryReader, demBEGasEThermoPhysicsH2TGD);

make2ChemistryReader(demPLEGasEThermoPhysicsH2TGD);
make2ChemistryReaderType(foam2ChemistryReader, demPLEGasEThermoPhysicsH2TGD);

make2ChemistryReader(demCEAGasEThermoPhysicsH2TGD);
make2ChemistryReaderType(foam2ChemistryReader, demCEAGasEThermoPhysicsH2TGD);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
