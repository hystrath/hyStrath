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

\*---------------------------------------------------------------------------*/

#include "patchLayer.H"
#include "transform.H"


namespace Foam
{
// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

patchLayer::patchLayer
(
    const polyMesh& mesh,
    const scalar& offset,
    const List<vector>& boundaryPoints, // 4 points
    const vector& normal,
    const vector& centroid
)
:
    mesh_(mesh),
    offset_(offset),
    boundaryPoints_(boundaryPoints),
    nF_(normal),
    centroid_(centroid)
{
    if(boundaryPoints_.size() != 4)
    {
        FatalErrorIn ("patchLayer.C")
            << nl << "BoundaryPoints has to be equal to 4 points. " << nl
            << abort(FatalError);
    }    
    
    createLayer();
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

patchLayer::~patchLayer()
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
void patchLayer::createLayer()
{
    createBox();
    findCells();
}

void patchLayer::findCells()
{
    DynamicList<label> cells(0);

    forAll(mesh_.cells(), c)
    {
        boundedBox bb = cellToBoundBox(c);
        
        if(bb_.justOverlaps(bb))
        {
            cells.append(c);
        }
    }
    
    cells.shrink();
    
    cells_.setSize(cells.size());
    cells_.transfer(cells);
}

void patchLayer::createBox()
{
    t1_ = (boundaryPoints_[0] - boundaryPoints_[1])/
            mag(boundaryPoints_[0] - boundaryPoints_[1]);
    t2_ = nF_ ^ t1_;
    t2_ /= Foam::mag(t2_);

    scalar xMax = 0.0;
    scalar xMin = GREAT;

    scalar yMax = 0.0;
    scalar yMin = GREAT;
    
    forAll(boundaryPoints_, i)
    {
        scalar x = (boundaryPoints_[i]-centroid_) & t1_;
        scalar y = (boundaryPoints_[i]-centroid_) & t2_;
        
        if(x > xMax)
        {
            xMax = x;
        }
        if(x < xMin)
        {
            xMin = x;
        }  
        
        if(y > yMax)
        {
            yMax = y;
        }
        if(y < yMin)
        {
            yMin = y;
        }         
    }
    
    vector min(xMin, yMin, 0);
    vector max(xMax, yMax, offset_);
    
    boundedBox bb (min, max);
    
    bb_ = bb;
}

void patchLayer::transformPointToNewCS(vector& point)
{
    vector newPoint =  vector
                        (
                            ((point - centroid_) & t1_),     
                            ((point - centroid_) & t2_),
                            ((point - centroid_) & nF_)
                        );
    
    point = newPoint;
}

boundedBox patchLayer::cellToBoundBox(const label& cellI)
{
    pointField points = boundingCellPoints(cellI);
    
    boundedBox bb(points, false);

    return bb;
}


pointField patchLayer::boundingCellPoints(const label& cellI)
{
    const labelList& points = mesh_.cellPoints()[cellI];
    const labelList& faces = mesh_.cells()[cellI];

    label sizeOfList = points.size() + faces.size();

    pointField vectorPoints(sizeOfList, vector::zero);

    label counter = 0;

    forAll(points, p)
    {
        vectorPoints[counter] = mesh_.points()[points[p]];
        counter++;
    }

    forAll(faces, f)
    {
        vectorPoints[counter] = mesh_.faceCentres()[faces[f]];
        counter++;
    }

    // transform
    forAll(vectorPoints, i)
    {
        transformPointToNewCS(vectorPoints[i]);
    }
    

    return vectorPoints;
}

boundedBox& patchLayer::bb()
{
    return bb_;
}

}

// ************************************************************************* //
