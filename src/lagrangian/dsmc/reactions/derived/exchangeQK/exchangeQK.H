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
    exchangeQK

Description


SourceFiles
    exchangeQK.C


\*---------------------------------------------------------------------------*/

#ifndef exchangeQK_H
#define exchangeQK_H

#include "dsmcReaction.H"
#include "vector.H"
#include "fvMesh.H"
#include "volFields.H"
#include "dsmcCloud.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class exchangeQK Declaration
\*---------------------------------------------------------------------------*/

class exchangeQK
:
   public virtual dsmcReaction
{

protected:

    // Protected data

        //- Sub-dictionary containing the reaction properties
        dictionary propsDict_;
        
        //- Position of the molecule in the list of reactants
        label posMolReactant_;
        
        //- Ids of the products
        labelList productIdsExchange_;
        
        //- String to print for the exchange reaction
        word exchangeStr_;
        
        label nTotExchangeReactions_;
        label nExchangeReactionsPerTimeStep_;
       
        scalar heatOfReactionExchangeJoules_;
        scalar aCoeff_, bCoeff_;

        //- Total volume of the numerical domain (NB: should not be there)
        scalar volume_;


    // Protected Member Functions

        virtual void setProperties();
        
        void testExchange
        (
            const dsmcParcel& p,
            const scalar translationalEnergy,
            const scalar omegaPQ,
            scalar& totalReactionProbability,
            scalar& reactionProbability
        );
        
        void exchange(dsmcParcel& p, dsmcParcel& q);


public:

    //- Runtime type information
    TypeName("exchangeQK");

    // Constructors

        //- Construct from components
        exchangeQK
        (
            Time& t,
            dsmcCloud& cloud,
            const dictionary& dict
        );
    
    
    // Destructor
    
        ~exchangeQK();
    
    
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

        virtual inline label nReactionsPerTimeStep() const
        {
            return nExchangeReactionsPerTimeStep_;
        }

        virtual void outputResults(const label& counterIndex);

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
