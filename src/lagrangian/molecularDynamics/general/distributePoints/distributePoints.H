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


Description

    distributing points in a boundedBox
    - Randomly - straightforward
    - Uniformly - the problem of distributing a number of points
            uniformly in a rectangle or cuboid is something which
            has not yet been determined properly in mathematics. 
            
    
    Note this: replaces the obsolete "randomPointsInSquareZone"


SourceFiles
    distributePointsI.H
    distributePoints.C
    distributePointsIO.C

\*---------------------------------------------------------------------------*/

#ifndef distributePoints_H
#define distributePoints_H

#include "IFstream.H"
#include "OFstream.H"
#include "polyMesh.H"
#include "cachedRandomMD.H"
#include "Time.H"
#include "boundedBox.H" 
#include "distribution.H" 

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class distributePoints Declaration
\*---------------------------------------------------------------------------*/

class distributePoints
:
    public boundedBox
{
    // Private data

        cachedRandomMD& rndGen_;

    // Private Member Functions
    
public:

    // Static data members

    // Constructors

    
     distributePoints
     (
         cachedRandomMD& rndGen
     );
    
     distributePoints
     (
         const boundedBox& bb,         
         cachedRandomMD& rndGen
     );     
     
//         //- Construct null
//         distributePoints
// 		(
// 			const label rndSeed,
// 			const label rndCacheSize
// 		);
// 
//         //- Construct
//         distributePoints
//         (
//             const boundedBox& bb,
// 			const label rndSeed,
// 			const label rndCacheSize
//         );

    // Selectors

    // Destructor

        ~distributePoints();

    // Member Functions

        void setBoundedBox(const boundedBox& bb);
        vector randomPoint();
		List<vector> uniform(const label& nPoints);
        void testUniform();
        void outputProperties(const List<vector>& r);
            
        // Access

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

#endif

// ************************************************************************* //
