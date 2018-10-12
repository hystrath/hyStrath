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

Description

\*---------------------------------------------------------------------------*/

#include "polyMesh.H"

namespace Foam
{


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


//- Null constructor
inline faceInfo::faceInfo()
:
//     List<label>(),
    face_(-1),
//     connectedPoints_(),
//     connectedCells_(),
    interfaceFace_(false),
    acceptedFace_(false),
    zoneFace_(false),
    boundaryFace_(false),
    cyclicFace_(false),
    procFace_(false)
{}


//- Construct from mesh and face
inline faceInfo::faceInfo
(
//     const polyMesh& mesh,
    const label& face
)
:
//     List<label>(mesh.faces()[face]),
    face_(face),
//     connectedPoints_(mesh.faces()[face]),
//     connectedCells_(),
    interfaceFace_(false),
    acceptedFace_(false),
    zoneFace_(false),
    boundaryFace_(false),
    cyclicFace_(false),
    procFace_(false)

{}


// //- Construct from components
// inline faceInfo::faceInfo
// (
// //     List<label>& points,
//     const label& face,
//     const labelList& points
// 
// )
// :
// //     List<label>(points),
//     face_(face),
//     connectedPoints_(points),
// //     connectedCells_(),
//     interfaceFace_(false),
//     acceptedFace_(false),
//     zoneFace_(false),
//     boundaryFace_(false),
//     cyclicFace_(false),
//     procFace_(false)
// {}



// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


inline void faceInfo::setFace(const label& face)
{
    face_ = face;
}

// inline void faceInfo::setConnectedPoints(const labelList& points)
// {
//     connectedPoints_ = points;
// }


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


inline void faceInfo::setFaceInfo
(
    const label& face
//     const labelList& points
)
{
    setFace(face);

//     setConnectedPoints(points);
}


// inline void faceInfo::setFaceInfo
// (
//     const polyMesh& mesh,
//     const label& face
// )
// {
//     face_ = face;
// 
//     connectedPoints_ = mesh.faces()[face];
// }




/*
inline void faceInfo::setFaceInfo
(
    PtrList<pointInfo>& pointData,
    const label& face
)
{
    face_ = face;

    transfer(pointData);
}*/


// inline void faceInfo::updateFaceInfo
// (
//     List<pointInfo>& pointData
// )
// {
//     checkAcceptance(pointData);
// }


// inline void faceInfo::setConnectedCells
// (
//     const polyMesh& mesh
// )
// {
//     DynamicList<label> connectedCells(0);
// 
//     const label& cellA = mesh.faceOwner()[face_];
// 
//     const label& cellB = mesh.faceNeighbour()[face_];
// 
//     connectedCells.append(cellA);
// 
//     if(cellB != -1)
//     {
//         connectedCells.append(cellB);
//     }
// 
//     connectedCells_ = connectedCells.shrink();
// }



// not used
inline void faceInfo::updateFaceBoundaryInfo
(
    const pointInfoList& pointData,
    const labelList& pointList
)
{
    if(!isBoundaryFace())
    {
        label noOfBoundaryPoints = 0;

        forAll(pointList, p)
        {
            if(pointData[pointList[p]].isBoundaryPoint())
            {
                noOfBoundaryPoints++;
            }
        }

        if(noOfBoundaryPoints == pointList.size())
        {
            faceOnBoundary();
        }
    }

    if(Pstream::parRun())
    {
        if(
            (!isProcBoundaryFace()) &&
            (isBoundaryFace())
          )
        {
            label noOfProcPoints = 0;

            forAll(pointList, p)
            {
                if(pointData[pointList[p]].isProcBoundaryPoint())
                {
                    noOfProcPoints++;
                }
            }

            if(noOfProcPoints == pointList.size())
            {
                faceOnProcBoundary();
            }
        }
    }

    if(
        (!isCyclicBoundaryFace()) &&
        (isBoundaryFace())
      )
    {
        label noOfCycPoints = 0;

        forAll(pointList, p)
        {
            if(pointData[pointList[p]].isCyclicBoundaryPoint())
            {
                noOfCycPoints++;
            }
        }

        if(noOfCycPoints == pointList.size())
        {
            faceOnCyclicBoundary();
        }
    }
}


// inline void faceInfo::checkAcceptance
// (
//     List<pointInfo>& pointData
// )
// {
//     if(!isAcceptedFace())
//     {
//         label noOfAcceptedPoints = 0;
// 
//         forAll(connectedPoints_, cP)
//         {
//             if(pointData[connectedPoints_[cP]].isAcceptedPoint())
//             {
//                 noOfAcceptedPoints++;
//             }
//         }
// 
//         if( noOfAcceptedPoints ==  connectedPoints_.size())
//         {
//             acceptFace();
// //             faceInZone();
//             faceOutOfZone();
//         }
//     }
// 
// 
// }

inline void faceInfo::faceOnInterface()
{
    interfaceFace_ = true;
}

inline void faceInfo::faceInZone()
{
    zoneFace_ = true;
}

inline void faceInfo::faceOutOfZone()
{
    zoneFace_ = false;
}

inline void faceInfo::acceptFace()
{
    acceptedFace_ = true;
}

inline void faceInfo::resetAcceptedFace()
{
    acceptedFace_ = false;
}

inline void faceInfo::faceOnBoundary()
{
    boundaryFace_ = true;
}

inline void faceInfo::faceOnProcBoundary()
{
    procFace_ = true;
}

inline void faceInfo::faceOnCyclicBoundary()
{
    cyclicFace_ = true;
}



// - Access

inline const label& faceInfo::faceLabel() const
{
    return face_;
}

// inline const labelList& faceInfo::connectedPoints() const
// {
//     return connectedPoints_;
// }

// inline const labelList& faceInfo::connectedCells(const polyMesh& mesh)
// {
//     setConnectedCells(mesh);
// 
//     return connectedCells_;
// }

inline const bool& faceInfo::isInterfaceFace() const
{
    return interfaceFace_;
}

inline const bool& faceInfo::isZoneFace() const
{
    return zoneFace_;
}

inline const bool& faceInfo::isAcceptedFace() const
{
    return acceptedFace_;
}

inline const bool& faceInfo::isBoundaryFace() const
{
    return boundaryFace_;
}

inline const bool& faceInfo::isProcBoundaryFace() const
{
    return procFace_;
}

inline const bool& faceInfo::isCyclicBoundaryFace() const
{
    return cyclicFace_;
}



// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //



} // End namespace Foam

// ************************************************************************* //
