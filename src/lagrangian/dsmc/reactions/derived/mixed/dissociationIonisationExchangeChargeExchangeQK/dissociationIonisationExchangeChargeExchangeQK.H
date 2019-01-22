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
    dissociationIonisationExchangeChargeExchangeQK

Description


SourceFiles
    dissociationIonisationExchangeChargeExchangeQK.C


\*---------------------------------------------------------------------------*/

#ifndef dissociationIonisationExchangeChargeExchangeQK_H
#define dissociationIonisationExchangeChargeExchangeQK_H

#include "dissociationQK.H"
#include "ionisationQK.H"
#include "exchangeQK.H"
#include "chargeExchangeQK.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
        Class dissociationIonisationExchangeChargeExchangeQK Declaration
\*---------------------------------------------------------------------------*/

class dissociationIonisationExchangeChargeExchangeQK
:
   public dissociationQK,
   public ionisationQK,
   public exchangeQK,
   public chargeExchangeQK
{

protected:

    // Protected data

        //- Sub-dictionary containing the reaction properties
        //dictionary propsDict_;
        

    // Protected Member Functions

        virtual void setProperties();


public:

    //- Runtime type information
    TypeName("dissociationIonisationExchangeChargeExchangeQK");

    // Constructors

        //- Construct from components
        dissociationIonisationExchangeChargeExchangeQK
        (
            Time& t,
            dsmcCloud& cloud,
            const dictionary& dict
        );
    
    
    // Destructor
    
        ~dissociationIonisationExchangeChargeExchangeQK();
    
    
    // Member Functions

        //- Initial configuration
        void initialConfiguration();

        //- Apply a reaction between parcel p and q
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

        virtual label nReactionsPerTimeStep() const;

        virtual void outputResults(const label& counterIndex);

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
