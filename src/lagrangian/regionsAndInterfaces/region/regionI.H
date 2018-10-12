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
    region

Description

\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#include "polyMesh.H"

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

inline void region::acceptPointForTraversing
(
    const label& traversePoint
)
{
    if(findIndex(pointsToTraverse_, traversePoint) == -1)
    {
        pointsToTraverse_.append(traversePoint);
    }
}


inline bool region::continueTraverse()
{
    pointsToTraverse_.shrink();

    if(pointsToTraverse_.size() == 0)
    {
        return false;
    }
    else
    {
        return true;
    }
}

inline bool region::contProcTraverse()
{
    if(contProcTraverse_)
    {
        return true;
    }
    else
    {
        return false;
    }
}


inline void region::traverseZone
(
    const polyMesh& mesh,
    pointInfoList& pointData,
    const label& originPoint,
    const label& originPointIndex
)
{
    traverse
    (
        mesh,
        pointData,
        originPoint,
        originPoint,
        originPointIndex
    );

    while(continueTraverse())
    {
        traversePoints_.clear();

        traversePoints_ = pointsToTraverse_;

        pointsToTraverse_.clear();

        forAll(traversePoints_, p)
        {
            traverse
            (
                mesh,
                pointData,
                originPoint,
                traversePoints_[p],
                originPointIndex
            );
        }
    }

    pointData.unVisitPoints();
}


inline void region::traverse
(
    const polyMesh& mesh,
    pointInfoList& pointData,
    const label& originPointLabel,
    const label& traverseNode,
    const label& originPointIndex
)
{
//     Info << " Traversing from point: " << startNode << endl;

    const labelList& neighbourPoints = mesh.pointPoints()[traverseNode];

//     Info << "  Neighbour points: " << neighbourPoints << endl;

    forAll(neighbourPoints, pointI)
    {
        const label& aP = neighbourPoints[pointI];

        if(pointData[aP].isZonePoint())
        {
//             if(!pointData[aP].isInterfacePoint())
//             {
                // - update point from start node

    //             Info << "   update point: " << aP << endl;

                const scalar& RsqrRef = pointData[traverseNode].Rsqr();

                const point& originPoint = pointData[traverseNode].originPoint();

                pointData[aP].updatePointInfo(originPoint, RsqrRef);

                const label& procNo = pointData[traverseNode].procNo();

                if(Pstream::parRun())
                {
                    pointData[aP].updateStartingProcNo(procNo);
                }

                if(!pointData[aP].isVisitedPoint())
                {
    //                 Info << "   visit point: " << aP << endl;

                    pointData[aP].visitPoint();

                    if(!pointData[aP].isAcceptedPoint())
                    {
    //                     Info << "   test point: " << aP << endl;
                        pointData[aP].testPoint();

                        if(pointData[aP].isAcceptedPoint())
                        {
                            acceptPointForTraversing(aP);

                            handleBoundaryPoint
                            (
                                mesh,
                                pointData, 
                                aP, // boundary point
                                originPoint,
                                originPointLabel,
                                originPointIndex,
                                procNo
                            );
                        }
                    }
                    else
                    {
                        acceptPointForTraversing(aP);

                        handleBoundaryPoint
                        (
                            mesh,
                            pointData, 
                            aP, // boundary point
                            originPoint,
                            originPointLabel,
                            originPointIndex,
                            procNo
                        );
                    }
                }
//             }
        }
    }
}

inline void region::handleBoundaryPoint
(
    const polyMesh& mesh,
    const pointInfoList& pointData,
    const label& bPoint,
    const point& originPoint,
    const label& originPointLabel,
    const label& originPointIndex,
    const label& procNo
)
{
    // - append boundary points
    if(Pstream::parRun())
    {
        if(
            (pointData[bPoint].isBoundaryPoint()) &&
            (!pointData[bPoint].isCyclicBoundaryPoint()) &&
            (!pointData[bPoint].isSwitchedProc()) &&
            (!pointData[bPoint].isInterfacePoint()) &&
            (pointData[bPoint].isProcBoundaryPoint())
          )
        {
            const labelList& patchList = pointData[bPoint].patchList();

            if(patchList.size() == 1)
            {
                const label& patchNo = patchList[0];

                if(findIndex(originPointpatches_[originPointIndex], patchNo) == -1)
                {
                    originPointpatches_[originPointIndex].append(patchNo);

                    const polyPatch& patch = mesh.boundaryMesh()[patchNo];
                    const labelList& meshPoints = patch.meshPoints();
                    const label& bPointIndex = findIndex(meshPoints, bPoint);

                    originPoints_[patchNo][bPointIndex].append(originPoint);
                    procNos_[patchNo][bPointIndex].append(procNo);
                }
            }
        }
    }

    if(
        (pointData[bPoint].isBoundaryPoint()) &&
        (!pointData[bPoint].isCyclicBoundaryPoint()) &&
        (!pointData[bPoint].isProcBoundaryPoint())
    )
    {
//         Info << "WARNING: External boundary has been reached" << endl;
// //                  << " Since the boundary is not a cyclic nor a processor boundary"
// //                  << " then the overlap Region will not fit within the domain \n"
// //                  << "   1. Reduce the overlap region size, or \n"
// //                  << "   2. Increase the domain size"
// //                  << endl;
    }
}


inline void region::handleProcBoundaries
(
    const polyMesh& mesh,
    pointInfoList& pointData
)
{

    // - visualisation of communicating points on processor boundaries 
//     DynamicList<label> procBoundaryPoints(0);
// 
//     forAll(originPoints_, patch)
//     {
//         if(originPoints_[patch].size() > 0)
//         {
//             forAll(originPoints_[patch], p)
//             {
//                 if(originPoints_[patch][p].size() > 0)
//                 {
//                     const polyPatch& patchI = mesh.boundaryMesh()[patch];
//                     const labelList& meshPoints = patchI.meshPoints();
// 
//                     procBoundaryPoints.append(meshPoints[p]);
//                 }
//             }
//         }
//     }
// 
//     procBoundaryPoints.shrink();
// 
//     pointsToVisualise_.append(procBoundaryPoints);


    Pout << "\t Sending data" << endl;

    //sending data to neighbouring processors
    forAll(mesh.boundaryMesh(), patchI)
    {
        const polyPatch& patch = mesh.boundaryMesh()[patchI];

        if (isA<processorPolyPatch>(patch))
        {
            // Note Rsqr remains the same for each subRegion.
            // Hence we don't need to pass Rsqrs across processors.

            List< List<point> > originPoints(originPoints_[patchI].size());
            List< List<label> > procNos(procNos_[patchI].size());

            forAll(originPoints_[patchI], p)
            {
                const DynamicList<point>& originPointList = originPoints_[patchI][p];
                const DynamicList<label>& procNosList = procNos_[patchI][p];

                originPoints[p].setSize(originPointList.size());
                procNos[p].setSize(procNosList.size());

                forAll(originPointList, oP)
                {
                    originPoints[p][oP] = originPointList[oP];
                    procNos[p][oP] = procNosList[oP];
                }

            }

            const processorPolyPatch& procPatch =
                    refCast<const processorPolyPatch>(patch);

            {
                OPstream toNeighbour(Pstream::commsTypes::blocking, procPatch.neighbProcNo());

                toNeighbour << originPoints << procNos;
            }
        }
    }

    Pout << "\t Receiving data" << endl;

    const label& patchListSize = mesh.boundaryMesh().size();

    List< List< List<point> > > procOriginPoints(patchListSize);
    List< List< List<label> > > procProcNos(patchListSize);

    //receiving data from neighbouring processors
    forAll(mesh.boundaryMesh(), patchI)
    {
        const polyPatch& patch = mesh.boundaryMesh()[patchI];

        if (isA<processorPolyPatch>(patch))
        {
            const processorPolyPatch& procPatch =
                refCast<const processorPolyPatch>(patch);

            List< List<point> > originPoints;
            List< List<label> > procNos;

            {
                IPstream fromNeighbour(Pstream::commsTypes::blocking, procPatch.neighbProcNo());

                fromNeighbour >> originPoints >> procNos;
            }

            procOriginPoints[patchI] = originPoints;
            procProcNos[patchI] = procNos;
        }
    }


    //clear lists
    forAll(originPoints_, oP)
    {
        forAll(originPoints_[oP], p)
        {
            originPoints_[oP][p].clear();
            procNos_[oP][p].clear();
        }
    }


    Pout << "\t Traversing from processor boundary" << endl;

    forAll(procOriginPoints, patch)
    {
        const List< List<point> >& originPatchPoints = procOriginPoints[patch];

        const List< List<label> >& procPatchNos = procProcNos[patch];

        const labelList& globalPatchPoints = correspPatchPoints_[patch];

        originPointpatches_.clear();

        originPointpatches_.setSize(globalPatchPoints.size());

        forAll(originPointpatches_, oP)
        {
            originPointpatches_[oP].setSize(0);
        }

        pointData.switchProcPoints(globalPatchPoints);

        forAll(originPatchPoints, p)
        {
            const List<point>& originPoints = originPatchPoints[p];

            const List<label>& procNos = procPatchNos[p];

            const label& patchPointI = globalPatchPoints[p];

            forAll(originPoints, oP)
            {
                const point& newOriginPoint = originPoints[oP];

                const label& newProcNo = procNos[oP];

                if(newProcNo != Pstream::myProcNo())
                {
                    pointData[patchPointI].updateStartingPointInfo(newOriginPoint, Rsqr_);

                    pointData[patchPointI].updateStartingProcNo(newProcNo);

                    traverseZone
                    (
                        mesh,
                        pointData,
                        patchPointI,
                        p
                    );
                }

            }
        }

        pointData.unSwitchProcPoints();
    }

    // shrink lists
    forAll(originPoints_, oP)
    {
        forAll(originPoints_[oP], p)
        {
            originPoints_[oP][p].shrink();
            procNos_[oP][p].shrink();
        }
    }
}



inline void region::correspPatchPoints
(
    const polyMesh& mesh,
    const polyPatch& patch,
    List<label>& patchPoints,
    const List<label>& patchPointFaces,
    const List<label>& facePointIndex
)
{
    forAll(patchPoints, p)
    {
        patchPoints[p] = correspondingPatchPoint
                         (
                            mesh,
                            patchPointFaces[p],
                            facePointIndex[p],
                            patch
                         );
    }
}


inline const label& region::correspondingPatchPoint
(
    const polyMesh& mesh,
    const label& faceIndex,
    const label& facePointIndex,
    const polyPatch& patch
) const
{

    const labelList& globalFaceI = mesh.faces()[patch.start() + faceIndex];

//     Pout << "global points: " << globalFaceI << endl;
//     const label& faceISize = faceI.size();
//     List<label> newPointList(faceISize);
//     newPointList[0] = faceI[0];
//     for (label i = 1; i <= faceISize; i++)
//     {
//        newPointList[i] = faceI[faceISize - i];
//     }

    if(facePointIndex == 0)
    {
        const label& correspPoint = globalFaceI[0];

        return correspPoint;
    }
    else
    {
        const label& faceSize = globalFaceI.size();

        const label& correspPoint = globalFaceI[faceSize - facePointIndex];

        return correspPoint;
    }

//     return newPointList[facePointIndex];
//     Pout << "oldPatchPoint: " << oldPatchPoint <<
//         "newPatchPoint: " << newPointList[facePointIndex] << endl;
}

/*
inline void pointInfoList::findGlobalPoint
(
    const polyMesh& mesh,
    const polyPatch& patch,
    const label& localPoint
)
{
    const labelListList& pointFacesList = patch.pointFaces();
    const faceList& localFaces = patch.localFaces();

    const label& faceI = pointFacesList[localPoint][0];

//     Pout << "faceI: " << faceI << endl;

    const label& pointFaceIndex = findIndex(localFaces[faceI], localPoint);

    const label globalFace = patch.start() + faceI;

//     Pout << "globalFace: " << globalFace << endl;

//     const labelList& pointList = mesh.faces()[globalFace];

//     const label& globalPoint = pointList[pointFaceIndex];

//     Pout << "Local Point " << localPoint 
//          << " Global Point: " << globalPoint << endl;


//     label startFaceI = 0;
// 
//     for(label i = 0; i < patch.size(); i++)
//     {
//         label patchFaceI = i + startFaceI;
// 
//         label meshFaceI = patch.start() + patchFaceI;
//     }

}*/






inline void region::checkAllProcs(const polyMesh& mesh)
{
    // This function checks whether the traverse across processor boundaries should terminate.
    // If one processor terminates before the others have finished, there will be errors of transmitting/receiving data.

//     Pout <<  "myProcNumber: " << Pstream::myProcNo() << endl;

    label pointsToTraverse = 0;

    forAll(originPoints_, oP)
    {
        forAll(originPoints_[oP], p)
        {
            pointsToTraverse += originPoints_[oP][p].size();
        }
    }

    if(procTraverse_[Pstream::myProcNo()])
    {
        if(pointsToTraverse == 0)
        {
            procTraverse_[Pstream::myProcNo()] = false;
        }
        else
        {
            procTraverse_[Pstream::myProcNo()] = true;
        }
    }

//     Pout <<  "procTraverse Before: " << procTraverse_ << endl;

//     Pout <<  "originPoints_ : " << originPoints_ << endl;

//     Pout << "sending check data: " << endl;

    // - sending data
    for (int i = 0; i < Pstream::nProcs(); i++)
    {
        if(i != Pstream::myProcNo())
        {
            const int proc = i;

            {
                OPstream toNeighbour(Pstream::commsTypes::blocking, proc);

                toNeighbour << procTraverse_;
            }
        }
    }

//     Pout << "receiving check data" << endl;

    // - receiving data
    for (int i = 0; i < Pstream::nProcs(); i++)
    {
        if(i != Pstream::myProcNo())
        {
            List<bool> procTraverse;

            const int proc = i;

            {
                IPstream fromNeighbour(Pstream::commsTypes::blocking, proc);

                fromNeighbour >> procTraverse;
            }

            forAll(procTraverse, p)
            {
                if(p != Pstream::myProcNo())
                {
                    if(!procTraverse[p])
                    {
                        procTraverse_[p] = false;
                    }
                }
            }
        }
    }

    Pout <<  "\t Processor traverse schedule: " << procTraverse_ << endl;

    label countProcs = 0;

    forAll(procTraverse_, p)
    {
        if(procTraverse_[p])
        {
            countProcs++;
        }
    }

    if(countProcs == 0)
    {
        contProcTraverse_ = false;
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// inline void region::setInterface(const labelList& faces)
// {
//     interfaceFaces_ = faces;
// }

// inline void region::setInterfaceName(const word& name)
// {
//     interfaceName_ = name;
// }


// - Access

inline const labelList& region::acceptedPoints() const
{
    return acceptedPoints_;
}

inline const labelList& region::regionCells() const
{
    return regionCells_;
}

inline const labelList& region::interfaceFaces() const
{
    return interfaceFaces_;
}

inline const labelList& region::interfacePoints() const
{
    return interfacePoints_;
}

inline const word& region::regionName() const
{
     return regionName_;
}

inline const word& region::interfaceName() const
{
    return interfaceName_;
}

// inline const DynamicList < DynamicList <label> >& region::pointsToVisualise() const
// {
//     return pointsToVisualise_;
// }

// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
