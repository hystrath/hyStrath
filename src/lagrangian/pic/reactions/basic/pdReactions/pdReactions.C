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

#include "pdReactions.H"
#include "pdCloud.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
//- Null Constructor
pdReactions::pdReactions
(
    Time& t,
    const polyMesh& mesh
)
:
    time_(t),
    chemReactDict_
    (
        IOobject
        (
            "chemReactDict",
            time_.system(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        )
    ),
    nReactions_(-1),
    reactionsList_(),
    reactions_(),
    pairAddressing_(),
    counter_(0)
{

}

//- Constructor for mdFOAM
pdReactions::pdReactions
(
    Time& t,
    const polyMesh& mesh,
    pdCloud& cloud
)
:
    time_(t),
    chemReactDict_
    (
        IOobject
        (
            "chemReactDict",
            time_.system(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    ),
    nReactions_(0),
    reactionsList_(chemReactDict_.lookup("reactions")),
    reactionNames_(reactionsList_.size()),
    reactionIds_(reactionsList_.size()),
    reactions_(reactionsList_.size()),
    counter_(0)
{

    Info << nl << "Creating pdReactions" << nl << endl;

    //- state pdReactions

    if(reactions_.size() > 0 )
    {
        forAll(reactions_, r)
        {
            const entry& pdReactionsI = reactionsList_[r];
            const dictionary& pdReactionsIDict = pdReactionsI.dict();

            reactions_[r] = autoPtr<pdReaction>
            (
                pdReaction::New(time_, cloud, pdReactionsIDict)
            );

            reactionNames_[r] = reactions_[r]->type();
            reactionIds_[r] = r;

            nReactions_++;
        }

        Info << "number of reactions created: " << nReactions_ << endl;
    }
    else
    {
        Info << "WARNING: there are no reactions." << endl;
    }

    pairAddressing_.setSize(cloud.typeIdList().size());

    forAll(pairAddressing_, p)
    {
        pairAddressing_[p].setSize(cloud.typeIdList().size(), -1);
    }
}

pdReactions::~pdReactions()
{
    //Info << "pdReactions Destructor" << endl;
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

//- initial configuration
//- call this function after the pdCloud is completely initialised
void pdReactions::initialConfiguration()
{
    forAll(reactions_, r)
    {
        reactions_[r]->initialConfiguration();
    }

    // set pair addressing

    forAll(pairAddressing_, i)
    {
        forAll(pairAddressing_[i], j)
        {
            label noOfReactionModelsPerPair = 0;

            forAll(reactions_, r)
            {
                if(reactions_[r]->tryReactMolecules(i, j))
                {
//                     Info << "r = " << r << endl;

                    pairAddressing_[i][j] = r;

//                     Info << "pairAddressing_[i][j] = " << pairAddressing_[i][j] << endl;

                    noOfReactionModelsPerPair++;
                }
            }

            if(noOfReactionModelsPerPair > 1)
            {
                FatalErrorIn("pdReactions::initialConfiguration()")
                    << "There is more than one reaction model specified for the typeId pair: "
                    << i << " and " << j
                    << exit(FatalError);
            }
        }
    }

    Info << "reactionNames: " << reactionNames_ << endl;

//     Info << "pair addressing: " << pairAddressing_ << endl;
}

//- output of data
//- required for output to screen any required information
void pdReactions::outputData()
{
    counter_++;

    forAll(reactions_, r)
    {
         reactions_[r]->outputResults(counter_);
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
