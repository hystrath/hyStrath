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
    associativeIonisationQK

Description


SourceFiles
    associativeIonisationQK.C


\*---------------------------------------------------------------------------*/

#ifndef associativeIonisationQK_H
#define associativeIonisationQK_H

#include "dsmcReaction.H"
#include "vector.H"
#include "fvMesh.H"
#include "volFields.H"
#include "dsmcCloud.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class associativeIonisationQK Declaration
\*---------------------------------------------------------------------------*/

class associativeIonisationQK
:
   public virtual dsmcReaction
{

protected:

    // Protected data

        //- Sub-dictionary containing the reaction properties
        dictionary propsDict_;
        
        //- Id of the intermediate molecule
        label intermediateMoleculeId_;
        
        //- Ids of the products
        labelList productIdsAssociativeIonisation_;
        
        //- Switch to determine whether or not the reaction is a forward or 
        //  reverse associative ionisation reaction
        bool forwardAssociativeIonisation_;
        
        //- Position of the ionised molecule in the reactant or product list
        label posIonisedMol_;
        
        //- String to print for the associative ionisation reaction
        word associativeIonisationStr_;
        
        label nTotAssociativeIonisationReactions_;
        label nAssociativeIonisationReactionsPerTimeStep_;
       
        scalarList heatOfReactionAssociativeIonisationJoules_;

        //- Total volume of the numerical domain (NB: should not be there)
        scalar volume_;


    // Protected Member Functions

        virtual void setProperties();
        
        void testForwardAssociativeIonisation
        (
            const dsmcParcel& p,
            const dsmcParcel& q,
            const scalar translationalEnergy,
            scalar& collisionEnergy,
            scalar& totalReactionProbability,
            scalar& reactionProbability
        );
        
        void testReverseAssociativeIonisation
        (
            const dsmcParcel& p,
            const dsmcParcel& q,
            const scalar translationalEnergy,
            scalar& collisionEnergy,
            scalar& totalReactionProbability,
            scalar& reactionProbability
        );
        
        void forwardAssociativeIonisation
        (
            dsmcParcel& p,
            dsmcParcel& q,
            scalar collisionEnergy
        );
        
        void reverseAssociativeIonisation
        (
            dsmcParcel& p,
            dsmcParcel& q,
            scalar collisionEnergy
        );


public:

    //- Runtime type information
    TypeName("associativeIonisationQK");

    // Constructors

        //- Construct from components
        associativeIonisationQK
        (
            Time& t,
            dsmcCloud& cloud,
            const dictionary& dict
        );
    
    
    // Destructor
    
        ~associativeIonisationQK();
    
    
    // Member Functions

        //- initial configuration - read in properties from chemReactDict
        void initialConfiguration();

        //- apply a reaction between parcel p and q
        void reaction
        (
            dsmcParcel& p,
            dsmcParcel& q,
            const DynamicList<label>& candidateList,
            const List<DynamicList<label>>& candidateSubList,
            const label& candidateP,
            const List<label>& whichSubCell
        );

        void reaction(dsmcParcel& p, dsmcParcel& q);

        //- Return true if both particles belong to this reaction model
        virtual bool tryReactMolecules
        (
            const label& typeIdP, 
            const label& typeIdQ
        );
        

      // Access

        virtual inline label nReactionsPerTimeStep() const
        {
            return nAssociativeIonisationReactionsPerTimeStep_;
        }

        virtual void outputResults(const label& counterIndex);

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
