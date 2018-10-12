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
    pointInfoList

Description

\*----------------------------------------------------------------------------*/

#include "pointInfoList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Null constructor
pointInfoList::pointInfoList()
:
    List<pointInfo> ()
{}

// constructor from list
pointInfoList::pointInfoList
(
    const List<pointInfo>& pointInfos 
)
:
    List<pointInfo> (pointInfos)
{}

// constructor from components
pointInfoList::pointInfoList
(
    const polyMesh& mesh
)
:
    List<pointInfo> (mesh.nPoints())
{
    setPoints(mesh);

    setBoundaryPoints(mesh);
}



// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

// Ostream& operator<<(Ostream& os, const pointInfoList& p)
// {
//     os << static_cast<const List<pointInfo>&>(p);
//     return os;
// }
// 
// 
// Istream& operator>>(Istream& is, pointInfoList& p)
// {
//     is >> static_cast<List<pointInfo>&>(p);
//     return is;
// }


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

pointInfoList::~pointInfoList()
{}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void pointInfoList::setBoundaryPoints( const polyMesh& mesh )
{
    forAll(mesh.boundaryMesh(), patchI)
    {
        //- points
        const polyPatch& patch = mesh.boundaryMesh()[patchI];

        const labelList& meshPoints = patch.meshPoints();

        forAll(meshPoints, i)
        {
            const label& globalPoint = meshPoints[i];

            operator[](globalPoint).pointOnBoundary();

            if (isA<cyclicPolyPatch>(patch))
            {
                operator[](globalPoint).pointOnCyclicBoundary();
            }

            if (isA<processorPolyPatch>(patch))
            {
                operator[](globalPoint).pointOnProcBoundary();
            }

            const labelList& faceList = mesh.pointFaces()[globalPoint];

//             Pout << "faceList: " << faceList << endl;

            DynamicList<label> patchIndices(0);

            forAll(faceList, f)
            {
                const label& faceI = faceList[f];

                // make sure this works for global faces rather than local faces
                const label& patchIndex = mesh.boundaryMesh().whichPatch(faceI); 

                if(patchIndex >= 0)
                {
                    if(findIndex(patchIndices, patchIndex) == -1)
                    {
                        patchIndices.append(patchIndex);
                    }
                }
            }

            operator[](globalPoint).setPatchList(patchIndices.shrink());
        }
    }
}



void pointInfoList::setPoints( const polyMesh& mesh )
{
    forAll(*this, p)
    {
        const point& pointI = mesh.points()[p];

        operator[](p).setPointInfo(pointI, label(p));
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void pointInfoList::clearInterfacePoints()
{
    forAll(*this, p)
    {
        operator[](p).resetInterfacePoint();
    }
}

void pointInfoList::unAcceptPoints()
{
    forAll(*this, p)
    {
        operator[](p).resetAcceptedPoint();
    }
}

void pointInfoList::clearZonePoints()
{
    forAll(*this, p)
    {
        operator[](p).pointOutOfZone();
    }
}

void pointInfoList::unVisitPoints()
{
    forAll(*this, p)
    {
        operator[](p).unVisitPoint();
    }
}

void pointInfoList::unSwitchProcPoints()
{
    forAll(*this, p)
    {
        operator[](p).unSwitchProc();
    }
}

void pointInfoList::switchProcPoints(const labelList& points)
{
    forAll(points, p)
    {
        operator[](points[p]).switchProc();
    }
}

// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
