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
    ionisationQK

Description


SourceFiles
    ionisationQK.C


\*---------------------------------------------------------------------------*/

#ifndef ionisationQK_H
#define ionisationQK_H

#include "dsmcReaction.H"
#include "vector.H"
#include "fvMesh.H"
#include "volFields.H"
#include "dsmcCloud.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class ionisationQK Declaration
\*---------------------------------------------------------------------------*/

class ionisationQK
:
   public virtual dsmcReaction
{

protected:

    // Protected data

        //- Sub-dictionary containing the reaction properties
        dictionary propsDict_;
        
        //- Ids of the products
        List<labelList> productIdsIonisation_;
        
        //- String to print for the ionisation reaction
        wordList ionisationStr_;
        
        labelList nTotIonisationReactions_;
        labelList nIonisationReactionsPerTimeStep_;
       
        scalarList heatOfReactionIonisationJoules_;

        //- Total volume of the numerical domain (NB: should not be there)
        scalar volume_;


    // Protected Member Functions

        virtual void setProperties();
        
        void testIonisation
        (
            const dsmcParcel& p,
            const scalar translationalEnergy,
            const label nReac,
            scalar& collisionEnergy,
            scalar& totalReactionProbability,
            scalar& reactionProbability
        );
        
        void ioniseParticleByPartner
        (
            dsmcParcel& p,
            dsmcParcel& q,
            const label nR,
            scalar collisionEnergy
        );


public:

    //- Runtime type information
    TypeName("ionisationQK");

    // Constructors

        //- Construct from components
        ionisationQK
        (
            Time& t,
            dsmcCloud& cloud,
            const dictionary& dict
        );
    
    
    // Destructor
    
        ~ionisationQK();
    
    
    // Member Functions

        //- initial configuration - read in properties from chemReactDict
        void initialConfiguration();

        //- apply a reaction between parcel p and q
        void reaction
        (
            dsmcParcel& p,
            dsmcParcel& q,
            const DynamicList<label>& candidateList,
            const List<DynamicList<label> >& candidateSubList,
            const label& candidateP,
            const List<label>& whichSubCell
        );

        void reaction(dsmcParcel& p, dsmcParcel& q);

        //- Return true if both particles belong to this reaction model
        virtual bool tryReactMolecules
        (
            const label &typeIdP, 
            const label& typeIdQ
        );
        

      // Access

        virtual label nReactionsPerTimeStep() const
        {
            return std::accumulate
            (
                nIonisationReactionsPerTimeStep_.begin(),
                nIonisationReactionsPerTimeStep_.end(),
                0
            );
        }

        virtual void outputResults(const label& counterIndex);

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
