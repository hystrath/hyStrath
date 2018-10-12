/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2007 OpenCFD Ltd.
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

Description

\*---------------------------------------------------------------------------*/

#include "pdReaction.H"
#include "IFstream.H"
#include "graph.H"
#include "pdCloud.H"


namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(pdReaction, 0);

defineRunTimeSelectionTable(pdReaction, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
pdReaction::pdReaction
(
    Time& t,
    pdCloud& cloud,
    const dictionary& dict
)
:
    mesh_(cloud.mesh()),
    cloud_(cloud),
    nTotReactions_(0),
    reactWithLists_(false)
{

}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

autoPtr<pdReaction> pdReaction::New
(
    Time& t,
    pdCloud& cloud,
    const dictionary& dict
)
{
    word pdReactionName
    (
        dict.lookup("reactionModel")
    );

    Info<< "Selecting the reaction model "
         << pdReactionName << endl;

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(pdReactionName);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalError
            << "pdReaction::New(const dictionary&) : " << endl
            << "    unknown pd reaction model type "
            << pdReactionName
            << ", constructor not in hash table" << endl << endl
            << "    Valid reaction types are :" << endl;
        Info<< dictionaryConstructorTablePtr_->toc() << abort(FatalError);
    }

    return autoPtr<pdReaction>
    (
        cstrIter()(t, cloud, dict)
    );
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

pdReaction::~pdReaction()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


const label& pdReaction::nTotReactions() const
{
    return nTotReactions_;
}

const label& pdReaction::nReactionsPerTimeStep() const
{
    return nReactionsPerTimeStep_;
}

label& pdReaction::nReactionsPerTimeStep()
{
    return nReactionsPerTimeStep_;
}
const pdCloud& pdReaction::cloud() const
{
    return cloud_;
}


const bool& pdReaction::reactWithLists() const
{
    return reactWithLists_;
}

} // End namespace Foam

// ************************************************************************* //
