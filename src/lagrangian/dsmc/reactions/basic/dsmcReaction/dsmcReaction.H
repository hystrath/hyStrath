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
    dsmcReaction

Description
    Basic/abstract class of a state controller


SourceFiles
    dsmcReactionI.H
    dsmcReaction.C
    dsmcReactionIO.C

\*---------------------------------------------------------------------------*/

#ifndef dsmcReaction_H
#define dsmcReaction_H

#include <numeric>

#include "IOdictionary.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "vector.H"
#include "volFields.H"
#include "Random.H"
#include "dsmcParcel.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class dsmcReaction Declaration
\*---------------------------------------------------------------------------*/

class dsmcReaction
{

protected:

    // Protected data

        const fvMesh& mesh_;

        dsmcCloud& cloud_;

        bool reactWithLists_;
        
        //- Name of the reaction
        word reactionName_;
        
        //- Ids of the reactants
        labelList reactantIds_;
        
        //- Types of the reactants
        labelList reactantTypes_;
        
        //- Switch to allow the chemical reaction to take place (default: true)
        //  Set to False to measure reaction rates in heat bath simulations 
        bool allowSplitting_;
        
        //- Switch to write the reaction rate in the log file (default: false)
        bool writeRatesToTerminal_;
        
        //- Switch to decide whether the collision is performed within the 
        //  reaction class or in the appropriate binary collision class
        bool relax_;
        

    // Private Member Functions
    
        //-
        virtual void setProperties();



public:

    //- Runtime type information
    TypeName("dsmcReaction");

    // Declare runtime constructor selection table
        declareRunTimeSelectionTable
        (
            autoPtr,
            dsmcReaction,
            dictionary,
            (
                Time& t,
                dsmcCloud& cloud,
                const dictionary& dict
            ),
            (t, cloud, dict)
        );
   
    // Constructors

        //- Construct from components
        dsmcReaction
        (
            Time& t,
            dsmcCloud& cloud,
            const dictionary& dict
        );


    // Selectors

        static autoPtr<dsmcReaction> New
        (
            Time& t,
            dsmcCloud& cloud,
            const dictionary& dict
        );


    // Destructor

        virtual ~dsmcReaction();

    // Member Functions

        virtual void initialConfiguration() = 0;

        virtual void reaction(dsmcParcel& p, dsmcParcel& q) = 0;

        virtual void reaction
        (
            dsmcParcel& p,
            dsmcParcel& q,
            const DynamicList<label>& candidateList,
            const List<DynamicList<label>>& candidateSubList,
            const label& candidateP,
            const List<label>& whichSubCell
        ) = 0;
        
        const bool& relax() const;

        virtual bool tryReactMolecules
        (
            const label& typeIdP,
            const label& typeIdQ
        ) = 0;

        virtual void outputResults(const label& counterIndex) = 0;
        
        labelList decreasing_sort_indices(const scalarList &v);

      // Access
      
        const dsmcCloud& cloud() const;

        virtual inline label nReactionsPerTimeStep() const = 0;

        const bool& reactWithLists() const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
