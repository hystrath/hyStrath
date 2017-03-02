/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
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

InClass
    Foam::make2Reaction

Description
    Macros for instantiating reactions on given thermo packages

\*---------------------------------------------------------------------------*/

#ifndef make2Reaction_H
#define make2Reaction_H

#include "Reaction2.H"

#include "Irreversible2Reaction.H"
#include "Reversible2Reaction.H"
#include "NonEquilibriumReversible2Reaction.H"

#include "multiThermo.H"

#include "constantTransport.H"
#include "SutherlandEuckenTransport.H"
#include "BlottnerEuckenTransport.H"
#include "powerLawEuckenTransport.H"
#include "CEATransport.H"

#include "janafThermo.H"
#include "decoupledEnergyModesThermo.H"
#include "perfectGas.H"

#include "polynomialTransport.H"
#include "hPolynomialThermo.H"
#include "icoPolynomial.H"

#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define makeReaction(MultiThermo, ReactionType, ReactionRate)                 \
                                                                              \
    typedef Reaction2<MultiThermo> Reaction2##MultiThermo;                    \
                                                                              \
    typedef ReactionType<Reaction2, MultiThermo, ReactionRate>                \
        ReactionType##MultiThermo##ReactionRate;                              \
                                                                              \
    template<>                                                                \
    const word ReactionType##MultiThermo##ReactionRate::typeName              \
    (                                                                         \
        ReactionType::typeName_()                                             \
      + ReactionRate::type()                                                  \
      + Reaction2##MultiThermo::typeName_()                                   \
    );                                                                        \
                                                                              \
    addToRunTimeSelectionTable                                                \
    (                                                                         \
        Reaction2##MultiThermo,                                               \
        ReactionType##MultiThermo##ReactionRate,                              \
        Istream                                                               \
    );                                                                        \
                                                                              \
    addToRunTimeSelectionTable                                                \
    (                                                                         \
        Reaction2##MultiThermo,                                               \
        ReactionType##MultiThermo##ReactionRate,                              \
        dictionary                                                            \
    );


#define makePressureDependentReaction\
(                                                                             \
    MultiThermo,                                                              \
    Reaction2,                                                                \
    PressureDependentReactionRate,                                            \
    ReactionRate,                                                             \
    FallOffFunction                                                           \
)                                                                             \
                                                                              \
    typedef PressureDependentReactionRate<ReactionRate, FallOffFunction>      \
        PressureDependentReactionRate##ReactionRate##FallOffFunction;         \
                                                                              \
    makeReaction                                                              \
    (                                                                         \
        MultiThermo,                                                          \
        Reaction2,                                                            \
        PressureDependentReactionRate##ReactionRate##FallOffFunction          \
    )


#define makeIRReactions(MultiThermo, ReactionRate)                            \
                                                                              \
    makeReaction(MultiThermo, Irreversible2Reaction, ReactionRate)            \
                                                                              \
    makeReaction(MultiThermo, Reversible2Reaction, ReactionRate)


#define makeIRNReactions(MultiThermo, ReactionRate)                           \
                                                                              \
    makeIRReactions(MultiThermo, ReactionRate)                                \
                                                                              \
    makeReaction(MultiThermo, NonEquilibriumReversible2Reaction, ReactionRate)


#define makePressureDependentReactions(MultiThermo, ReactionRate, FallOffFunction) \
                                                                              \
    makePressureDependentReaction                                             \
    (                                                                         \
        MultiThermo,                                                          \
        Irreversible2Reaction,                                                \
        FallOffReactionRate,                                                  \
        ReactionRate,                                                         \
        FallOffFunction                                                       \
    )                                                                         \
                                                                              \
    makePressureDependentReaction                                             \
    (                                                                         \
        MultiThermo,                                                          \
        Reversible2Reaction,                                                  \
        FallOffReactionRate,                                                  \
        ReactionRate,                                                         \
        FallOffFunction                                                       \
    )                                                                         \
                                                                              \
    makePressureDependentReaction                                             \
    (                                                                         \
        MultiThermo,                                                          \
        Irreversible2Reaction,                                                \
        ChemicallyActivatedReactionRate,                                      \
        ReactionRate,                                                         \
        FallOffFunction                                                       \
    )                                                                         \
                                                                              \
    makePressureDependentReaction                                             \
    (                                                                         \
        MultiThermo,                                                          \
        Reversible2Reaction,                                                  \
        ChemicallyActivatedReactionRate,                                      \
        ReactionRate,                                                         \
        FallOffFunction                                                       \
    )


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
