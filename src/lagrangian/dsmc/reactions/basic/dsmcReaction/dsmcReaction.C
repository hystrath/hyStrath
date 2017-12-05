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

#include "dsmcReaction.H"
#include "IFstream.H"
#include "graph.H"
#include "dsmcCloud.H"


namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(dsmcReaction, 0);

defineRunTimeSelectionTable(dsmcReaction, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
dsmcReaction::dsmcReaction
(
    Time& t,
    dsmcCloud& cloud,
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

autoPtr<dsmcReaction> dsmcReaction::New
(
    Time& t,
    dsmcCloud& cloud,
    const dictionary& dict
)
{
    word dsmcReactionName
    (
        dict.lookup("reactionModel")
    );

    Info<< "Selecting the reaction model "
         << dsmcReactionName << endl;

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(dsmcReactionName);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalError
            << "dsmcReaction::New(const dictionary&) : " << endl
            << "    unknown dsmc reaction model type "
            << dsmcReactionName
            << ", constructor not in hash table" << endl << endl
            << "    Valid reaction types are :" << endl;
        Info<< dictionaryConstructorTablePtr_->toc() << abort(FatalError);
    }

    return autoPtr<dsmcReaction>
    (
        cstrIter()(t, cloud, dict)
    );
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

dsmcReaction::~dsmcReaction()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


const label& dsmcReaction::nTotReactions() const
{
    return nTotReactions_;
}

const label& dsmcReaction::nReactionsPerTimeStep() const
{
    return nReactionsPerTimeStep_;
}

label& dsmcReaction::nReactionsPerTimeStep()
{
    return nReactionsPerTimeStep_;
}
const dsmcCloud& dsmcReaction::cloud() const
{
    return cloud_;
}


const bool& dsmcReaction::reactWithLists() const
{
    return reactWithLists_;
}

} // End namespace Foam

// ************************************************************************* //
