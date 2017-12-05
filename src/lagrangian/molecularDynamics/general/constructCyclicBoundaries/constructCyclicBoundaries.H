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
    // IMPORTANT RUNS IN SERIAL
    
    This class is used by the createCyclicBoundaries utility to generate information
    on the cyclic boundaries. This information is stored appropriately in 
    the case/system/ directory in a dictionary called "cyclicBoundaries" which is
    created in this class. This dictionary is then re-read when running the MD simulations
    in serial or parallel, The information allows the proper initialisation 
    of the "cyclicBoundaries" object. 

SourceFiles
    constructCyclicBoundariesI.H
    constructCyclicBoundaries.C
    constructCyclicBoundariesIO.C

\*---------------------------------------------------------------------------*/

#ifndef constructCyclicBoundaries_H
#define constructCyclicBoundaries_H

#include "IFstream.H"
#include "OFstream.H"
#include "polyMesh.H"
#include "Time.H"
#include "IOdictionary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class constructCyclicBoundaries Declaration
\*---------------------------------------------------------------------------*/

class constructCyclicBoundaries
{
    // Private data

        const polyMesh& mesh_;

        label nCyclicBoundaries_;

        List< List<vector> > boundaryPoints_;
        
        List<word> names_;
        
        List<label> nFaces_;

        List<vector> normal_;

        List<vector> centroid_;

        List<word> neighbourNames_;        
        
        
    // Private Member Functions

            void setPatchPoints
            (
                const label& patchId,
                const label& c
            );
            
            void setCyclicBoundaries();    


public:

    // Static data members

    // Constructors

        //- Construct 
        constructCyclicBoundaries
        (
            const polyMesh& mesh
        );
    

    // Selectors

    // Destructor

        ~constructCyclicBoundaries();

    // Member Functions

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
