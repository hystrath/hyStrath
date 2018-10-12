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
inline cellInform::cellInform()
:
//     List<faceInfo>(),
    cell_(-1),
//     connectedFaces_(),
    acceptedCell_(false),
    prefZoneCell_(false),
    zoneCell_(false),
    boundaryCell_(false)
{}


// //- Construct from mesh and cell
// inline cellInform::cellInform
// (
// //     const polyMesh& mesh,
//     const label& cell
// )
// :
// //     List<faceInfo>(mesh.cells()[cell].size()),
//     cell_(cell),
// //     connectedFaces_(mesh.cells()[cell]),
//     acceptedCell_(false),
//     prefZoneCell_(false),
//     zoneCell_(false),
//     boundaryCell_(false)
// {}


//- Construct from components
inline cellInform::cellInform
(
    const label& cell
//     const labelList& faces
)
:
//     List<faceInfo>(faceData),
    cell_(cell),
//     connectedFaces_(faces),
    acceptedCell_(false),
    prefZoneCell_(false),
    zoneCell_(false),
    boundaryCell_(false)
{}





// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

inline void cellInform::setCell(const label& cell)
{
    cell_ = cell;
}

// inline void cellInform::setConnectedFaces(const labelList& faces)
// {
//     connectedFaces_ = faces;
// }

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


inline void cellInform::setCellInfo
(
    const label& cell
//     const labelList& faces
)
{
    setCell(cell);

//     setConnectedFaces(faces);
}

/*
inline void cellInform::setCellInfo
(
    const polyMesh& mesh,
    const label& cell
)
{
    cell_ = cell;

    connectedFaces_ =mesh.cells()[cell];
}*/



// inline void cellInform::updateCellInfo()
// {
//     checkAcceptance();
// }


// inline void cellInform::checkAcceptance
// (
//     const List<faceInfo>& faceData,
//     const labelList& faces
// )
// {
//     if(!isAcceptedCell())
//     {
//         label noOfAcceptedFaces = 0;
// 
//         forAll(connectedFaces_, cF)
//         {
// 
//             if(faceData[connectedFaces_[cF]].isAcceptedFace())
//             {
//                 noOfAcceptedFaces++;
//             }
//         }
// 
//         if( noOfAcceptedFaces ==  connectedFaces_.size())
//         {
// //             cellOnBoundary();
//         }
//     }
// }

inline void cellInform::checkCellBoundary
(
    const faceInfoList& faceData,
    const labelList& faces
)
{
    forAll(faces, f)
    {
        if(!isBoundaryCell())
        {
            if(faceData[faces[f]].isBoundaryFace())
            {
                cellOnBoundary();
            }
        }
    }
}


// inline void cellInform::checkBoundary
// (
//     const faceInfoList& faceData
// )
// {
//     if(!isBoundaryCell())
//     {
//         forAll(connectedFaces_, cF)
//         {
//             if(faceData[connectedFaces_[cF]].isBoundaryFace())
//             {
//                 cellOnBoundary();
// 
//                 break;
//             }
//         }
//     }
// }


inline void cellInform::cellInZone()
{
    zoneCell_ = true;
}

inline void cellInform::cellOutOfZone()
{
    zoneCell_ = false;
}

inline void cellInform::cellInPrefZone()
{
    prefZoneCell_ = true;
}


inline void cellInform::cellOutOfPrefZone()
{
    prefZoneCell_ = false;
}


inline void cellInform::acceptCell()
{
    acceptedCell_ = true;
}

inline void cellInform::resetAcceptedCell()
{
    acceptedCell_ = false;
}

inline void cellInform::cellOnBoundary()
{
    boundaryCell_ = true;
}


// inline void cellInform::pointOnBoundary()
// {
//     boundaryPoint_ = true;
// }
// 
// inline void cellInform::pointOnProcBoundary()
// {
//     procBoundaryPoint_ = true;
// }
// 
// inline void cellInform::pointOnCyclicBoundary()
// {
//     cyclicBoundaryPoint_ = true;
// }
// 
// inline void cellInform::pointOnInterface()
// {
//     interfacePoint_ = true;
// }
// 
// inline void cellInform::resetInterfacePoint()
// {
//     interfacePoint_ = false;
// }



// inline void cellInform::switchDomain
// (
//     const polyMesh mesh,
//     const label& p
// )
// {
//     setCurrentPoint(p);
// 
//     switchedDomain_ = true;
// 
//     noOfSwitchedDomains_++;
// }


// - Access

inline const label& cellInform::cellLabel() const
{
    return cell_;
}

// inline const labelList& cellInform::connectedFaces() const
// {
//     return connectedFaces_;
// }

inline const bool& cellInform::isZoneCell() const
{
    return zoneCell_;
}

inline const bool& cellInform::isPrefZoneCell() const
{
    return prefZoneCell_;
}

inline const bool& cellInform::isAcceptedCell() const
{
    return acceptedCell_;
}

inline const bool& cellInform::isBoundaryCell() const
{
    return boundaryCell_;
}






// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //



} // End namespace Foam

// ************************************************************************* //
