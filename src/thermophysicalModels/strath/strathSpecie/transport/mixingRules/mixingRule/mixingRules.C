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

#include "makeMixingRule.H"

#include "molarMR.H"
//#include "GuptaMR.H"
#include "WilkeMR.H"
#include "ArmalySuttonMR.H"

#include "thermoPhysics2Types.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/* * * * * * * * * * * * * * * private static data * * * * * * * * * * * * * */

// Based on sensible enthalpy

makeMixingRule(molarMR, demConstGasHThermoPhysicsH2TGD);
//makeMixingRule(GuptaMR, demConstGasHThermoPhysicsH2TGD);
makeMixingRule(WilkeMR, demConstGasHThermoPhysicsH2TGD);
makeMixingRule(ArmalySuttonMR, demConstGasHThermoPhysicsH2TGD);

makeMixingRule(molarMR, demGasHThermoPhysicsH2TGD);
//makeMixingRule(GuptaMR, demGasHThermoPhysicsH2TGD);
makeMixingRule(WilkeMR, demGasHThermoPhysicsH2TGD);
makeMixingRule(ArmalySuttonMR, demGasHThermoPhysicsH2TGD);

makeMixingRule(molarMR, demBEGasHThermoPhysicsH2TGD);
//akeMixingRule(GuptaMR, demBEGasHThermoPhysicsH2TGD);
makeMixingRule(WilkeMR, demBEGasHThermoPhysicsH2TGD);
makeMixingRule(ArmalySuttonMR, demBEGasHThermoPhysicsH2TGD);

makeMixingRule(molarMR, demPLEGasHThermoPhysicsH2TGD);
//makeMixingRule(GuptaMR, demPLEGasHThermoPhysicsH2TGD);
makeMixingRule(WilkeMR, demPLEGasHThermoPhysicsH2TGD);
makeMixingRule(ArmalySuttonMR, demPLEGasHThermoPhysicsH2TGD);

makeMixingRule(molarMR, demCEAGasHThermoPhysicsH2TGD);
//makeMixingRule(GuptaMR, demCEAGasHThermoPhysicsH2TGD);
makeMixingRule(WilkeMR, demCEAGasHThermoPhysicsH2TGD);
makeMixingRule(ArmalySuttonMR, demCEAGasHThermoPhysicsH2TGD);


// Based on sensible internal energy

makeMixingRule(molarMR, demConstGasEThermoPhysicsH2TGD);
//makeMixingRule(GuptaMR, demConstGasEThermoPhysicsH2TGD);
makeMixingRule(WilkeMR, demConstGasEThermoPhysicsH2TGD);
makeMixingRule(ArmalySuttonMR, demConstGasEThermoPhysicsH2TGD);

makeMixingRule(molarMR, demGasEThermoPhysicsH2TGD);
//makeMixingRule(GuptaMR, demGasEThermoPhysicsH2TGD);
makeMixingRule(WilkeMR, demGasEThermoPhysicsH2TGD);
makeMixingRule(ArmalySuttonMR, demGasEThermoPhysicsH2TGD);

makeMixingRule(molarMR, demBEGasEThermoPhysicsH2TGD);
//makeMixingRule(GuptaMR, demBEGasEThermoPhysicsH2TGD);
makeMixingRule(WilkeMR, demBEGasEThermoPhysicsH2TGD);
makeMixingRule(ArmalySuttonMR, demBEGasEThermoPhysicsH2TGD);

makeMixingRule(molarMR, demPLEGasEThermoPhysicsH2TGD);
//makeMixingRule(GuptaMR, demPLEGasEThermoPhysicsH2TGD);
makeMixingRule(WilkeMR, demPLEGasEThermoPhysicsH2TGD);
makeMixingRule(ArmalySuttonMR, demPLEGasEThermoPhysicsH2TGD);

makeMixingRule(molarMR, demCEAGasEThermoPhysicsH2TGD);
//makeMixingRule(GuptaMR, demCEAGasEThermoPhysicsH2TGD);
makeMixingRule(WilkeMR, demCEAGasEThermoPhysicsH2TGD);
makeMixingRule(ArmalySuttonMR, demCEAGasEThermoPhysicsH2TGD);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
