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

    Generating a random point in a square zone on the mesh 

    as well in parallel

SourceFiles
    randomPointsInSquareZoneI.H
    randomPointsInSquareZone.C
    randomPointsInSquareZoneIO.C

\*---------------------------------------------------------------------------*/

#ifndef randomPointsInSquareZone_H
#define randomPointsInSquareZone_H

#include "IFstream.H"
#include "OFstream.H"
#include "polyMesh.H"
#include "Random.H"
#include "Time.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class randomPointsInSquareZone Declaration
\*---------------------------------------------------------------------------*/

class randomPointsInSquareZone
{
    // Private data

//         const polyMesh& mesh_;

//         label regionId_;

        scalar deltaX_;
        scalar deltaY_;
        scalar deltaZ_;
        scalar XMin_;
        scalar YMin_;
        scalar ZMin_;

        Random rndGen_;

    // Private Member Functions

    
public:

    // Static data members

    // Constructors

        //- Construct null
        randomPointsInSquareZone();

        //- Construct from xs and ys
        randomPointsInSquareZone
        (
            const polyMesh& mesh,
            const label& regionId
        );
    
        randomPointsInSquareZone
        (
            const boundBox& bb
        );

    // Selectors

    // Destructor

        ~randomPointsInSquareZone();

    // Member Functions

        void setInitialConfiguration
        (
            const boundBox& bb
        );

        void setInitialConfiguration
        (
            const polyMesh& mesh,
            const label& regionId
        );

        // Access

            vector randomPoint();



        // Check

        // Edit

        // Write


    // Member Operators

//         void operator=(const randomPointsInSquareZone&);

    // Friend Functions

    // Friend Operators

    // IOstream Operators

//         friend Ostream& operator<<(Ostream&, const randomPointsInSquareZone&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
