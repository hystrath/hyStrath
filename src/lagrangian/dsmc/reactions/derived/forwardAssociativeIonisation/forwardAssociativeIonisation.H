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
    forwardAssociativeIonisation

Description


SourceFiles
    forwardAssociativeIonisation.C


\*---------------------------------------------------------------------------*/

#ifndef forwardAssociativeIonisation_H
#define forwardAssociativeIonisation_H

#include "dsmcReaction.H"
#include "vector.H"
#include "fvMesh.H"
#include "volFields.H"
#include "dsmcCloud.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class forwardAssociativeIonisation Declaration
\*---------------------------------------------------------------------------*/

class forwardAssociativeIonisation
:
   public dsmcReaction
{

private:
    // Private data

        //- subDictionary containing the properties
        dictionary propsDict_;

        List<label> reactantIds_;
        
        label intermediateId_;

        List<label> associativeIonisationProductIds_;
        List<label> ionisationProductIds_;
        
        word reactionName_;
        
        scalar heatOfReactionRecombination_;
        scalar heatOfReactionIntermediateIonisation_;
        scalar heatOfReactionIonisation_;
        
        label nTotalIonisationReactions_;
        label nTotalAssociativeIonisationReactions_;
        
        label nIonisationReactionsPerTimeStep_;
        label nAssociativeIonisationReactionsPerTimeStep_;

        bool relax_;

        // remove if not used for heat bath
        bool allowSplitting_;
        bool writeRatesToTerminal_;
        label nReactions_;
        scalar volume_;
        List<scalar> numberDensities_;


    // Private Member Functions

        void setProperties();

public:

    //- Runtime type information
    TypeName("forwardAssociativeIonisation");

    // Constructors

        //- Construct from components
        forwardAssociativeIonisation
        (
	        Time& t,
            dsmcCloud& cloud,
            const dictionary& dict
        );
    
        // Destructor
        ~forwardAssociativeIonisation();
    
    // Member Functions

        // initial configuration - read in properties from chemReactDict
        void initialConfiguration();

        // apply a reaction between parcel p and q
        void reaction
        (
            dsmcParcel& p,
            dsmcParcel& q,
            const DynamicList<label>& candidateList,
            const List<DynamicList<label> >& candidateSubList,
            const label& candidateP,
            const List<label>& whichSubCell
        );

        void reaction
        (
            dsmcParcel& p,
            dsmcParcel& q
        );

        // return bool if molecules belong to this reaction model
        bool tryReactMolecules(const label& typeIdP, const label& typeIdQ);

        // Access

        // return bool for relax = true (run normal collision model)
        // relax = false (apply reaction to molecules)
        const bool& relax() const;

        void outputResults(const label& counterIndex);

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
