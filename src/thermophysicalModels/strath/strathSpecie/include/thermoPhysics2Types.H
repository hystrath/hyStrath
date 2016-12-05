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

Typedefs
    Foam::thermoPhysics2Types

Description
    Type definitions for thermo-physics models, multi-temperature models

\*---------------------------------------------------------------------------*/

#ifndef thermoPhysics2Types_H
#define thermoPhysics2Types_H

#include "advancedSpecie.H"
#include "perfect2Gas.H"
#include "decoupledEnergyModesThermo.H"

#include "multiThermo.H"
#include "sensible2Enthalpy.H"
#include "sensible2InternalEnergy.H"

#include "constantTransport.H"
#include "SutherlandEuckenTransport.H"
#include "BlottnerEuckenTransport.H"
#include "powerLawEuckenTransport.H"
#include "CEATransport.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    // thermo physics types based on sensibleEnthalpy
    
    typedef
    constantTransport
    <
        species::multiThermo
        <
            decoupledEnergyModesThermo
            <
                perfect2Gas<advancedSpecie>
            >,
            sensible2Enthalpy
        >
    > demConstGasHThermoPhysicsH2TGD;
    
    typedef
    SutherlandEuckenTransport
    <
        species::multiThermo
        <
            decoupledEnergyModesThermo
            <
                perfect2Gas<advancedSpecie>
            >,
            sensible2Enthalpy
        >
    > demGasHThermoPhysicsH2TGD;
    
    typedef
    BlottnerEuckenTransport
    <
        species::multiThermo
        <
            decoupledEnergyModesThermo
            <
                perfect2Gas<advancedSpecie>
            >,
            sensible2Enthalpy
        >
    > demBEGasHThermoPhysicsH2TGD;
    
    typedef
    powerLawEuckenTransport
    <
        species::multiThermo
        <
            decoupledEnergyModesThermo
            <
                perfect2Gas<advancedSpecie>
            >,
            sensible2Enthalpy
        >
    > demPLEGasHThermoPhysicsH2TGD;
    
    typedef
    CEATransport
    <
        species::multiThermo
        <
            decoupledEnergyModesThermo
            <
                perfect2Gas<advancedSpecie>
            >,
            sensible2Enthalpy
        >
    > demCEAGasHThermoPhysicsH2TGD;

    // thermo physics types based on sensibleInternalEnergy
    
    typedef
    constantTransport
    <
        species::multiThermo
        <
            decoupledEnergyModesThermo
            <
                perfect2Gas<advancedSpecie>
            >,
            sensible2InternalEnergy
        >
    > demConstGasEThermoPhysicsH2TGD;
    
    typedef
    SutherlandEuckenTransport
    <
        species::multiThermo
        <
            decoupledEnergyModesThermo
            <
                perfect2Gas<advancedSpecie>
            >,
            sensible2InternalEnergy
        >
    > demGasEThermoPhysicsH2TGD;
    
    typedef
    BlottnerEuckenTransport
    <
        species::multiThermo
        <
            decoupledEnergyModesThermo
            <
                perfect2Gas<advancedSpecie>
            >,
            sensible2InternalEnergy
        >
    > demBEGasEThermoPhysicsH2TGD;
    
    typedef
    powerLawEuckenTransport
    <
        species::multiThermo
        <
            decoupledEnergyModesThermo
            <
                perfect2Gas<advancedSpecie>
            >,
            sensible2InternalEnergy
        >
    > demPLEGasEThermoPhysicsH2TGD;
    
    typedef
    CEATransport
    <
        species::multiThermo
        <
            decoupledEnergyModesThermo
            <
                perfect2Gas<advancedSpecie>
            >,
            sensible2InternalEnergy
        >
    > demCEAGasEThermoPhysicsH2TGD;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
