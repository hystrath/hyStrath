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


// * * * * * * * * * * *  Protected Member functions * * * * * * * * * * * * //

void dsmcReaction::setProperties()
{
    if (cloud_.coordSystem().type() != "dsmcCartesian" && writeRatesToTerminal_)
    {
        WarningIn("dsmcReaction::setProperties()")
            << "For reaction named " << reactionName_ << nl
            << "Non Cartesian coordinate system, cannot print "
            << "reaction rates to the terminal." << nl << endl;
            writeRatesToTerminal_ = false;
    }
    
    forAll(reactantIds_, r)
    {
        reactantTypes_[r] = cloud_.constProps(reactantIds_[r]).type();
    }
}


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
    reactWithLists_(false),
    reactionName_(dict.name().substr(1)),
    reactantIds_(),
    reactantTypes_(),
    allowSplitting_
    (
        dict.lookupOrDefault<Switch>("allowSplitting", true)
    ),
    writeRatesToTerminal_
    (
        dict.lookupOrDefault<Switch>("writeRatesToTerminal", false)
    ),
    relax_(true)
{
    //- Reading in reactants
    const List<word> reactants(dict.lookup("reactants"));
    reactantIds_.setSize(reactants.size(), -1);
    reactantTypes_.setSize(reactantIds_.size(), -1);
    
    forAll(reactantIds_, r)
    {
        reactantIds_[r] = findIndex(cloud_.typeIdList(), reactants[r]);
        
        //- Check that reactants belong to the typeIdList as defined in 
        //  constant/dsmcProperties
        if (reactantIds_[r] == -1)
        {
            FatalErrorIn("dsmcReaction::setProperties()")
                << "For reaction named " << reactionName_ << nl
                << "Cannot find type id: " << reactants[r] << nl 
                << exit(FatalError);
        }
    }
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

const dsmcCloud& dsmcReaction::cloud() const
{
    return cloud_;
}


const bool& dsmcReaction::reactWithLists() const
{
    return reactWithLists_;
}


const bool& dsmcReaction::relax() const
{
    return relax_;
}


labelList dsmcReaction::decreasing_sort_indices(const scalarList &v)
{
    //- initialize original index locations
    labelList idx(v.size());
    std::iota(idx.begin(), idx.end(), 0);

    scalarList rnd(v.size(), 0.0);
    forAll(rnd, i)
    {
        rnd[i] = cloud_.rndGen().sample01<scalar>();
    }

    //- sort indices in decreasing order based on comparing values in v
    //  random shuffle for similar values
    std::sort
    (
        idx.begin(),
        idx.end(),
        [&v, &rnd](label i1, label i2)
        {
            if (v[i1] == v[i2])
            {
                return rnd[i1] > rnd[i2];
            }
            return v[i1] > v[i2];
        }
    );

    return idx;
}

} // End namespace Foam

// ************************************************************************* //
