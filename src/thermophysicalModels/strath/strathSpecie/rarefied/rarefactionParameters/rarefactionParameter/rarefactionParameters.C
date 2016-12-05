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

#include "makeRarefactionParameter.H"

#include "rarefied.H"

#include "thermoPhysics2Types.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/* * * * * * * * * * * * * * * private static data * * * * * * * * * * * * * */

// Based on sensible enthalpy

makeRarefactionParameter(rarefied, demConstGasHThermoPhysicsH2TGD);

makeRarefactionParameter(rarefied, demGasHThermoPhysicsH2TGD);

makeRarefactionParameter(rarefied, demBEGasHThermoPhysicsH2TGD);

makeRarefactionParameter(rarefied, demPLEGasHThermoPhysicsH2TGD);

makeRarefactionParameter(rarefied, demCEAGasHThermoPhysicsH2TGD);

// Based on sensible internal energy

makeRarefactionParameter(rarefied, demConstGasEThermoPhysicsH2TGD);

makeRarefactionParameter(rarefied, demGasEThermoPhysicsH2TGD);

makeRarefactionParameter(rarefied, demBEGasEThermoPhysicsH2TGD);

makeRarefactionParameter(rarefied, demPLEGasEThermoPhysicsH2TGD);

makeRarefactionParameter(rarefied, demCEAGasEThermoPhysicsH2TGD);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
