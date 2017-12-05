/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2005 OpenCFD Ltd.
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
    dsmcDynamicLoadBalancing

Description

SourceFiles
    dsmcDynamicLoadBalancingI.H
    dsmcDynamicLoadBalancing.C

\*---------------------------------------------------------------------------*/

#ifndef dsmcDynamicLoadBalancing_H
#define dsmcDynamicLoadBalancing_H

#include "vector.H"
#include "fvMesh.H"
#include "dsmcParcel.H"
#include "decompositionMethod.H"
#include "volFields.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class dsmcDynamicLoadBalancing Declaration
\*---------------------------------------------------------------------------*/
// class dsmcCloud;


class dsmcDynamicLoadBalancing
{
    // Private data

        Time& time_;
        
        const fvMesh& mesh_;

        dsmcCloud& cloud_;
        
        IOdictionary dsmcLoadBalanceDict_;
        
        //IOdictionary controlDict_;
        
        bool performBalance_;
        
        bool enableBalancing_;
        
        scalar originalEndTime_;
        
        scalar maxImbalance_;
        
//         label nProcs_;

    // Private Member Functions



public:

    friend class dsmcParcel;

    // Constructors

         //- Constructor
        dsmcDynamicLoadBalancing
        (
            Time& t,
            const polyMesh& mesh,
            dsmcCloud& cloud
        );

    // Destructor

        ~dsmcDynamicLoadBalancing();

    // Member Functions

        void update();
        
        void perform
        (
            const int noRefinement = 0
        );
        
        void copyPolyMeshToLatestTimeFolder() const;

        
        void updateProperties(/*const IOdictionary&*/);

        // Access
        
          inline bool performBalance() const
          {
              return performBalance_;
          }


        // Check

        // Edit

        // Write

    // Member Operators

    // Friend Functions

    // Friend Operators


    // IOstream Operators

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//#include "dsmcDynamicLoadBalancingI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
