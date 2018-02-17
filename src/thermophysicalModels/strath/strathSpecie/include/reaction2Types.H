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
    Foam::reaction2Types

Description
    Type definitions for reactions, multi-temperature models

\*---------------------------------------------------------------------------*/

#ifndef reaction2Types_H
#define reaction2Types_H

#include "thermoPhysics2Types.H"
#include "Reaction2.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    // sensible enthalpy based reactions
    
    typedef Reaction2<demConstGasHThermoPhysicsH2TGD> demConstGasHReactionH2TGD;
    
    typedef Reaction2<demGasHThermoPhysicsH2TGD> demGasHReactionH2TGD;
    
    typedef Reaction2<demBEGasHThermoPhysicsH2TGD> demBEGasHReactionH2TGD;
    
    typedef Reaction2<demPLEGasHThermoPhysicsH2TGD> demPLEGasHReactionH2TGD;
    
    typedef Reaction2<demCEAGasHThermoPhysicsH2TGD> demCEAGasHReactionH2TGD;

    // internal energy based reactions
    
    typedef Reaction2<demConstGasEThermoPhysicsH2TGD> demConstGasEReactionH2TGD;
    
    typedef Reaction2<demGasEThermoPhysicsH2TGD> demGasEReactionH2TGD;
    
    typedef Reaction2<demBEGasEThermoPhysicsH2TGD> demBEGasEReactionH2TGD;
    
    typedef Reaction2<demPLEGasEThermoPhysicsH2TGD> demPLEGasEReactionH2TGD;
    
    typedef Reaction2<demCEAGasEThermoPhysicsH2TGD> demCEAGasEReactionH2TGD;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
