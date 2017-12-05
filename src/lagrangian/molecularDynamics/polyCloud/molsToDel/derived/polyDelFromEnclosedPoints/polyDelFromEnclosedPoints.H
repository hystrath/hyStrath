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
    polyDelFromEnclosedPoints

Description
    Apply 2D points 

SourceFiles
    polyDelFromEnclosedPointsI.H
    polyDelFromEnclosedPoints.C
    polyDelFromEnclosedPointsIO.C

\*---------------------------------------------------------------------------*/

#ifndef polyDelFromEnclosedPoints_H
#define polyDelFromEnclosedPoints_H

#include "polyMolsToDeleteModel.H"
#include "vector.H"
#include "polyMesh.H"
#include "volFields.H"
#include "polyMoleculeCloud.H"
#include "boundBox.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class polyDelFromEnclosedPoints Declaration
\*---------------------------------------------------------------------------*/

class polyDelFromEnclosedPoints
:
	public polyMolsToDeleteModel
{

private:

    // Private data

        //- subDictionary containing the properties
        dictionary propsDict_;
        
        boundBox box_;
        
        List<label> molIds_;

        bool deleteExternalMols_;
        
        List<vector> points_;
        
        vector normal_;
        
        vector centre_;
        
        bool rotate_;
        
        // radians
        scalar theta_;
        
        
        word option_;
        
        
        List<vector> centreLine_;
        
        // outer list = size of points_;
        // inner list = 3 neighbours
        // inner list = centreLine_;
        List<List<vector> > c_;
        
        
        
        // local search
//         bool local_;
//         scalar radius_;
        
        
    // Private Member Functions

        void findMolsToDel();
        
        bool isPointWithinRegion(const vector& rI);
        
        bool isPointWithinTriangle
        (
            const vector& centre,
            const vector& r1,
            const vector& r2,
            const vector& rI
        );
        
        
public:

    //- Runtime type information
    TypeName("polyDelFromEnclosedPoints");

    // Constructors

        //- Construct from components
        polyDelFromEnclosedPoints
		(
            polyMoleculeCloud& molCloud,
            const dictionary& dict
		);

    // Destructor

 		~polyDelFromEnclosedPoints();


    // Member Functions

 
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
