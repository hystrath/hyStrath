/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012-2013 OpenFOAM Foundation
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

#include "reaction2Types.H"
#include "make2Reaction.H"

#include "Arrhenius2ReactionRate.H"
#include "subArrhenius2ReactionRate.H" // NEW VINCENT 08/11/16
#include "ImpactIonisationArrhenius2ReactionRate.H" // NEW VINCENT 03/06/16
#include "infiniteReactionRate.H"
#include "LandauTellerReactionRate.H"
#include "thirdBodyArrheniusReactionRate.H"

#include "ChemicallyActivatedReactionRate.H"
#include "JanevReactionRate.H"
#include "powerSeriesReactionRate.H"

#include "FallOffReactionRate.H"
#include "LindemannFallOffFunction.H"
#include "SRIFallOffFunction.H"
#include "TroeFallOffFunction.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define makeReactions(MultiThermo, Reaction)                                   \
                                                                               \
    defineTemplateTypeNameAndDebug(Reaction, 0);                               \
    defineTemplateRunTimeSelectionTable(Reaction, Istream);                    \
    defineTemplateRunTimeSelectionTable(Reaction, dictionary);                 \
                                                                               \
    makeIRNReactions(MultiThermo, Arrhenius2ReactionRate)                      \
    makeIRNReactions(MultiThermo, subArrhenius2ReactionRate)                   \
    makeIRNReactions(MultiThermo, ImpactIonisationArrhenius2ReactionRate)      \
    makeIRNReactions(MultiThermo, infiniteReactionRate)                        \
    makeIRNReactions(MultiThermo, LandauTellerReactionRate)                    \
    makeIRNReactions(MultiThermo, thirdBodyArrheniusReactionRate)              \
                                                                               \
    makeIRReactions(MultiThermo, JanevReactionRate)                            \
    makeIRReactions(MultiThermo, powerSeriesReactionRate)                      \
                                                                               \
    makePressureDependentReactions                                             \
    (                                                                          \
       MultiThermo,                                                            \
       Arrhenius2ReactionRate,                                                 \
       LindemannFallOffFunction                                                \
    )                                                                          \
                                                                               \
    makePressureDependentReactions                                             \
    (                                                                          \
       MultiThermo,                                                            \
       subArrhenius2ReactionRate,                                              \
       LindemannFallOffFunction                                                \
    )                                                                          \
                                                                               \
    makePressureDependentReactions                                             \
    (                                                                          \
       MultiThermo,                                                            \
       Arrhenius2ReactionRate,                                                 \
       TroeFallOffFunction                                                     \
    )                                                                          \
                                                                               \
    makePressureDependentReactions                                             \
    (                                                                          \
       MultiThermo,                                                            \
       Arrhenius2ReactionRate,                                                 \
       SRIFallOffFunction                                                      \
    )


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    // sensible enthalpy based reactions
    makeReactions(demConstGasHThermoPhysicsH2TGD, demConstGasHReactionH2TGD);
    
    makeReactions(demGasHThermoPhysicsH2TGD, demGasHReactionH2TGD);
    
    makeReactions(demBEGasHThermoPhysicsH2TGD, demBEGasHReactionH2TGD);
    
    makeReactions(demPLEGasHThermoPhysicsH2TGD, demPLEGasHReactionH2TGD);
    
    makeReactions(demCEAGasHThermoPhysicsH2TGD, demCEAGasHReactionH2TGD);
    
    // sensible internal energy based reactions
    makeReactions(demConstGasEThermoPhysicsH2TGD, demConstGasEReactionH2TGD);
    
    makeReactions(demGasEThermoPhysicsH2TGD, demGasEReactionH2TGD);
    
    makeReactions(demBEGasEThermoPhysicsH2TGD, demBEGasEReactionH2TGD);
    
    makeReactions(demPLEGasEThermoPhysicsH2TGD, demPLEGasEReactionH2TGD);
    
    makeReactions(demCEAGasEThermoPhysicsH2TGD, demCEAGasEReactionH2TGD);
}

// ************************************************************************* //
