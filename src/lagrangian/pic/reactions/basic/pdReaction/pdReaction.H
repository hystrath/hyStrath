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

Class
    pdReaction

Description
    Basic/abstract class of a state controller


SourceFiles
    pdReactionI.H
    pdReaction.C
    pdReactionIO.C

\*---------------------------------------------------------------------------*/

#ifndef pdReaction_H
#define pdReaction_H

#include "IOdictionary.H"
// #include "Time.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "vector.H"
#include "volFields.H"
#include "Random.H"
#include "pdParcel.H"
// #include "timeData.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class pdReaction Declaration
\*---------------------------------------------------------------------------*/

class pdReaction
{

protected:

    // Protected data


        const fvMesh& mesh_;

        pdCloud& cloud_;

        label nTotReactions_;

        label nReactionsPerTimeStep_;

        bool reactWithLists_;

      // Private Member Functions


public:

    //- Runtime type information
    TypeName("pdReaction");

    // Declare runtime constructor selection table
        declareRunTimeSelectionTable
        (
            autoPtr,
            pdReaction,
            dictionary,
            (
                Time& t,
                pdCloud& cloud,
                const dictionary& dict
            ),
            (t, cloud, dict)
        );
   
    // Constructors

        //- Construct from components
        pdReaction
        (
			Time& t,
            pdCloud& cloud,
            const dictionary& dict
        );


    // Selectors

        static autoPtr<pdReaction> New
        (
			Time& t,
            pdCloud& cloud,
            const dictionary& dict
        );


    // Destructor

        virtual ~pdReaction();

    // Member Functions

        virtual void initialConfiguration() = 0;

        virtual void reaction
        (
            pdParcel& p,
            pdParcel& q
        ) = 0;

        virtual void reaction
        (
            pdParcel& p,
            pdParcel& q,
            const DynamicList<label>& candidateList,
            const List<DynamicList<label> >& candidateSubList,
            const label& candidateP,
            const List<label>& whichSubCell
        ) = 0;

        virtual const bool& relax() const = 0;

        virtual bool tryReactMolecules(const label& typeIdP, const label& typeIdQ) = 0;

        virtual void outputResults(const label& counterIndex) = 0;

        // Access
        const pdCloud& cloud() const;

        const label& nTotReactions() const;

        const label& nReactionsPerTimeStep() const;
        label& nReactionsPerTimeStep();

        const bool& reactWithLists() const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
