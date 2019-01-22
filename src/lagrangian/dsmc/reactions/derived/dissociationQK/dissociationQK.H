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
    dissociationQK

Description


SourceFiles
    dissociationQK.C


\*---------------------------------------------------------------------------*/

#ifndef dissociationQK_H
#define dissociationQK_H

#include "dsmcReaction.H"
#include "vector.H"
#include "fvMesh.H"
#include "volFields.H"
#include "dsmcCloud.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class dissociationQK Declaration
\*---------------------------------------------------------------------------*/

class dissociationQK
:
   public virtual dsmcReaction
{

protected:

    // Protected data

        //- Sub-dictionary containing the reaction properties
        dictionary propsDict_;
        
        //- Ids of the products
        List<labelList> productIdsDissociation_;
        
        //- String to print for the dissociation reaction
        wordList dissociationStr_;
        
        labelList nTotDissociationReactions_;
        labelList nDissociationReactionsPerTimeStep_;
       
        List<scalarList> heatOfReactionDissociationJoules_;

        //- Total volume of the numerical domain (NB: should not be there)
        scalar volume_;


    // Protected Member Functions

        virtual void setProperties();
        
        void testDissociation
        (
            const dsmcParcel& p,
            const scalar translationalEnergy,
            label& vibModeDisso,
            scalar& collisionEnergy,
            scalar& totalReactionProbability,
            scalar& reactionProbability
        );
        
        void dissociateParticleByPartner
        (
            dsmcParcel& p,
            dsmcParcel& q,
            const label nR,
            const label vibModeDisso,
            scalar collisionEnergy
        );


public:

    //- Runtime type information
    TypeName("dissociationQK");

    // Constructors

        //- Construct from components
        dissociationQK
        (
            Time& t,
            dsmcCloud& cloud,
            const dictionary& dict
        );
    
    
    // Destructor
    
        ~dissociationQK();
    
    
    // Member Functions

        //- Initial configuration - read in properties from chemReactDict
        void initialConfiguration();

        //- Apply a reaction between parcel p and q
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
            const label& typeIdP,
            const label& typeIdQ
        );
        

      // Access

        virtual label nReactionsPerTimeStep() const
        {
            return std::accumulate
            (
                nDissociationReactionsPerTimeStep_.begin(),
                nDissociationReactionsPerTimeStep_.end(),
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
