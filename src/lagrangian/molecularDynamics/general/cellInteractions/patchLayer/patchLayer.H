/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::patchLayer

Description
    Creates a layer of cells defined by an offset distance from a patch
    
    Assumes patch is a flat plane and rectangular
    
    Patch is therefore defined by 4 points 

    Builds a boundedBox that starts from points on patch with an offset thickness 
    
    Any cells just overlapping with the boundedBox are selected

    [note] - still need to test this for inclined planes
    
SourceFiles
    patchLayerI.H
    patchLayer.C

\*---------------------------------------------------------------------------*/

#ifndef patchLayer_H
#define patchLayer_H


#include "polyMesh.H"
#include "vector.H"
#include "clock.H"
#include "boundedBox.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class patchLayer Declaration
\*---------------------------------------------------------------------------*/

class patchLayer
{
    // Private data
        const polyMesh& mesh_;
        const scalar& offset_;
        const List<vector>& boundaryPoints_;
        vector nF_;
        vector centroid_;
        labelList cells_;
        vector t1_;
        vector t2_;
        boundedBox bb_;
        
        void createLayer();
        void createBox();
        void findCells();
        void transformPointToNewCS(vector& point);
        boundedBox cellToBoundBox(const label& cellI);
        pointField boundingCellPoints(const label& cellI);
        
public:

    // Constructors

        //- Construct from components
        patchLayer
        (
            const polyMesh& mesh,
            const scalar& offset,
            const List<vector>& boundaryPoints,
            const vector& normal,
            const vector& centroid
        );

    //- Destructor
    ~patchLayer();


    // Member Functions

        // Access

		//- Contains point? (inside only)
        inline const labelList& cells() const;
        boundedBox& bb();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "patchLayerI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
