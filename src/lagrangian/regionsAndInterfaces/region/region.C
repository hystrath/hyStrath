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

\*----------------------------------------------------------------------------*/

#include "region.H"

// /home/graham/OpenFOAM/graham-1.4.1/run/gnemdFOAM/stanfordMixer_T_2.5_equalPotentials

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

Ostream& operator<<(Ostream& os, const region& sL)
{
    if (os.format() == IOstream::ASCII)
    {
        return os
            << sL.regionCells();
    }
    else
    {
        return os
            << sL.regionCells();
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// - construct null

region::region()
:
    Rsqr_(GREAT),
    interfacePoints_(),
    pointsToTraverse_(0),
    traversePoints_(0),
    regionCells_(),
    interfaceFaces_(),
    originPoints_(),
    procNos_(),
    originPointpatches_(),
    correspPatchPoints_(),
//     pointsToVisualise_(),
    acceptedPoints_(),
    contProcTraverse_(false),
    procTraverse_()

{}


// - construct from R

region::region
(
    const scalar& R
)
:
    Rsqr_(sqr(R)),
    interfacePoints_(),
    pointsToTraverse_(0),
    traversePoints_(0),
    regionCells_(),
    interfaceFaces_(),
    originPoints_(),
    procNos_(),
    originPointpatches_(),
    correspPatchPoints_(),
//     pointsToVisualise_(),
    acceptedPoints_(),
    contProcTraverse_(false),
    procTraverse_()
{}

// - construct from Rsqr and interfacePoints.

region::region
(
    const scalar& Rsqr,
    const labelList& interfacePoints
)
:
    Rsqr_(Rsqr),
    interfacePoints_(interfacePoints),
    pointsToTraverse_(0),
    traversePoints_(0),
    regionCells_(),
    interfaceFaces_(),
    originPoints_(),
    procNos_(),
    originPointpatches_(),
    correspPatchPoints_(),
//     pointsToVisualise_(),
    acceptedPoints_(),
    contProcTraverse_(false),
    procTraverse_()
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

region::~region()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void region::setRegion
(
    const polyMesh& mesh,
    const labelList& cellZone,
    const labelList& interface,
    const scalar& Rsqr,
    const word& name,
    pointInfoList& pointData,
    faceInfoList& faceData,
    cellInfoList& cellData
)
{
    //make cellsInZone
    cellData.cellsInZone(cellZone);

    setZonePoints(mesh, cellZone, pointData);

    //determine the corresponding patch points across processor boundaries
    setCorrPatchPoints(mesh, pointData);

    setRsqr(Rsqr);

    setInterface(mesh, interface, pointData);

    setRegionName(name);

    setRegionPoints(mesh, pointData);

    setRegionCells(mesh, pointData, cellData);

    setNewInterface(mesh, faceData, cellData);

    //reset data
    cellData.cellsOutOfZone(cellZone);
    pointData.clearInterfacePoints();
    pointData.unAcceptPoints();
    pointData.clearZonePoints();

}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void region::setZonePoints
(
    const polyMesh& mesh,
    const labelList& cellZone,
    pointInfoList& pointData
)
{
    // points in zone
    forAll(cellZone, cell)
    {
        const label& cellI = cellZone[cell];
        const labelList& pointList = mesh.cellPoints()[cellI];

        forAll(pointList, p)
        {
            pointData[pointList[p]].pointInZone();
        }
    }
}


void region::setRsqr(const scalar& Rsqr)
{
    Rsqr_ = Rsqr;
}


void region::setInterface
(
    const polyMesh& mesh,
    const labelList& interface,
    pointInfoList& pointData
)
{
    //assign interface points
    forAll(interface, face)
    {
        const labelList& pointList = mesh.faces()[interface[face]];

        forAll(pointList, p)
        {
            pointData[pointList[p]].pointOnInterface();
        }
    }

    //synchronise interface points which are on a processor boundaryMesh
    synchIntProcPoints(mesh, pointData);


    DynamicList<label> interfacePoints(0);

    forAll(pointData, p)
    {
        if( pointData[p].isInterfacePoint() )
        {
           interfacePoints.append(pointData[p].currentPointLabel());
        }
    }

    interfacePoints_ = interfacePoints.shrink();
}

//          Alternative to "originPointpatches_":

//     if(Pstream::parRun())
//     {
//
//         processorTraverse_.setSize(interfacePoints_);
// 
//        const label& procNo = Pstream::myProcNo();
// 
//         forAll(processorTraverse_, p)
//         {
//             processorTraverse_[p].setSize(Pstream::nProcs());
// 
//             forAll(processorTraverse_[p], procI)
//             {
//                 if(procI == procNo)
//                 {
//                     processorTraverse_[p][procI] = true;
//                 }
//                 else
//                 {
//                     processorTraverse_[p][procI] = false;
//                 }
//             }
//         }
//     }





// The most imp function of the region class:
// Traverse from each point on the interface in order to 
// accept all neighbouring points that lie within R.
void region::setRegionPoints
(
    const polyMesh& mesh,
    pointInfoList& pointData
)
{
    if(Pstream::parRun())
    {
        originPointpatches_.setSize(interfacePoints_.size());

        forAll(originPointpatches_, oP)
        {
            originPointpatches_[oP].setSize(0);
        }
    }

    // Determine accepted points internal to the domain
    forAll(interfacePoints_, intPoint)
    {
        const label& interfacePoint = interfacePoints_[intPoint];

        const point& iPoint = mesh.points()[interfacePoint];

        pointData[interfacePoint].updateStartingPointInfo(iPoint, Rsqr_);

        if(Pstream::parRun())
        {
            pointData[interfacePoint].updateStartingProcNo(Pstream::myProcNo());
        }

        // traverse from each of the interface points
        traverseZone(mesh, pointData, interfacePoint, intPoint);
    }

    handleBoundaries(mesh, pointData);

    setAcceptedPoints(pointData);
}


void region::handleBoundaries
(
    const polyMesh& mesh,
    pointInfoList& pointData
)
{
    if(Pstream::parRun())
    {
        contProcTraverse_ = true;

        procTraverse_.setSize(Pstream::nProcs());

        forAll(procTraverse_, p)
        {
            procTraverse_[p] = true;
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

        label terminationCounter = 0;

        // handle processor boundaries
        while(contProcTraverse())
        {
            handleProcBoundaries
            (
                mesh,
                pointData
            );

            checkAllProcs(mesh);

            terminationCounter++;

            if(terminationCounter >= Pstream::nProcs()*3)
            {
                contProcTraverse_ = false;

                Pout << "WARNING: looping terminated. Check domain-decomposition"
                     << endl;
            }
        }
//         pointsToVisualise_.shrink();
    }
}


void region::setRegionCells
(
    const polyMesh& mesh,
    pointInfoList& pointData,
    cellInfoList& cellData
)
{
    DynamicList<label> regionCells(0);

    forAll(pointData, p)
    {
        if(pointData[p].isAcceptedPoint())
        {
            const label& pointI = pointData[p].currentPointLabel();

            const labelList& cellList = mesh.pointCells()[pointI];

            forAll(cellList, c)
            {
                const label& cellI = cellList[c];

                if(cellData[cellI].isZoneCell())
                {
                    if(!cellData[cellI].isAcceptedCell())
                    {
                        cellData[cellI].acceptCell();

                        regionCells.append(cellI);
                    }
                }
            }
        }
    }

    regionCells_ = regionCells.shrink();

}

void region::setNewInterface
(
    const polyMesh& mesh,
    faceInfoList& faceData,
    cellInfoList& cellData
)
{
    DynamicList<label> remainingCells(0);

    forAll(cellData, c)
    {
        if(
            (!cellData[c].isAcceptedCell()) &&
            (cellData[c].isZoneCell())
          )
        {
            remainingCells.append(cellData[c].cellLabel());
        }
    }

    cellData.resetAcceptedCells(regionCells_);

    interface newInterface
    (
        mesh,
        faceData,
        cellData,
        remainingCells.shrink(),
        regionName_,
        regionCells_,
        "subLayerInterface"
    );

    interfaceFaces_ = newInterface.interfaceFaces();

    interfaceName_ = newInterface.interfaceName();
}

void region::setAcceptedPoints(pointInfoList& pointData)
{
    DynamicList<label> acceptedPoints(0);

    forAll(pointData, p)
    {
        if(pointData[p].isAcceptedPoint())
        {
            acceptedPoints.append(pointData[p].currentPointLabel());
        }
    }

    acceptedPoints_ = acceptedPoints.shrink();
}




// sets a list of corresponding points across processor boundaries 
// (this is done once before building the region)
void region::setCorrPatchPoints
(
    const polyMesh& mesh,
    const pointInfoList& pointData
)
{
    if(Pstream::parRun())
    {
        const label& patchListSize = mesh.boundaryMesh().size();

        correspPatchPoints_.setSize(patchListSize);

        //- sending
        forAll(mesh.boundaryMesh(), patchI)
        {
            const polyPatch& patch = mesh.boundaryMesh()[patchI];

            const labelList& meshPoints = patch.meshPoints();

            if (isA<processorPolyPatch>(patch))
            {
                const labelListList& pointFacesList = patch.pointFaces();
                const faceList& localFaces = patch.localFaces();

                DynamicList<label> patchPoints(0); // list of local points on the patch
                DynamicList<label> patchPointFaces(0); //list of face labels linked to the local patchPoints
                DynamicList<label> facePointIndex(0); // index of local point within faceList

                correspPatchPoints_[patchI].setSize(meshPoints.size());

                forAll(meshPoints, i)
                {
                    const label& localMeshPoint = patch.whichPoint(meshPoints[i]);
                    patchPoints.append(localMeshPoint);

                    const label& firstFace = pointFacesList[i][0];
                    patchPointFaces.append(firstFace); // 0 takes the first face on the list

                    facePointIndex.append
                    (
                        findIndex
                        (
                            localFaces[firstFace],
                            localMeshPoint
                        )
                    );
                }

                patchPoints.shrink();
                patchPointFaces.shrink();
                facePointIndex.shrink();

                const processorPolyPatch& procPatch =
                    refCast<const processorPolyPatch>(patch);

                {
                    OPstream toNeighbour(Pstream::commsTypes::blocking, procPatch.neighbProcNo());

                    toNeighbour << patchPoints << patchPointFaces << facePointIndex;
                }
            }
        }

        //-receiving
        forAll(mesh.boundaryMesh(), patchI)
        {
            const polyPatch& patch = mesh.boundaryMesh()[patchI];

            if (isA<processorPolyPatch>(patch))
            {
                const processorPolyPatch& procPatch =
                    refCast<const processorPolyPatch>(patch);

                List<label> patchPoints;
                List<label> patchPointFaces;
                List<label> facePointIndex;

                {
                    IPstream fromNeighbour(Pstream::commsTypes::blocking, procPatch.neighbProcNo());

                    fromNeighbour >> patchPoints >> patchPointFaces >> facePointIndex;
                }

                correspPatchPoints
                (
                    mesh,
                    patch,
                    patchPoints,
                    patchPointFaces,
                    facePointIndex
                );

                correspPatchPoints_[patchI] = patchPoints;
            }
        }


        originPoints_.setSize(patchListSize);
        procNos_.setSize(patchListSize);

        forAll(mesh.boundaryMesh(), patchI)
        {
            const polyPatch& patch = mesh.boundaryMesh()[patchI];
            const labelList& meshPoints = patch.meshPoints();

            originPoints_[patchI].setSize(meshPoints.size());
            procNos_[patchI].setSize(meshPoints.size());

            forAll(originPoints_[patchI], oP)
            {
                originPoints_[patchI][oP].setSize(0);
                procNos_[patchI][oP].setSize(0);
            }
        }
    }
}


void region::synchIntProcPoints
(
    const polyMesh& mesh,
    pointInfoList& pointData
)
{
    if(Pstream::parRun())
    {
        //- sending
        forAll(mesh.boundaryMesh(), patchI)
        {
            const polyPatch& patch = mesh.boundaryMesh()[patchI];

            if (isA<processorPolyPatch>(patch))
            {
                const labelList& meshPoints = patch.meshPoints();

                List<bool> intProcPoints(meshPoints.size());

                forAll(meshPoints, p)
                {
                    if(
                        (pointData[meshPoints[p]].isInterfacePoint()) &&
                        (pointData[meshPoints[p]].isProcBoundaryPoint())
                    )
                    {
                        intProcPoints[p] = true;
                    }
                    else
                    {
                        intProcPoints[p] = false;
                    }
                }

                const processorPolyPatch& procPatch =
                        refCast<const processorPolyPatch>(patch);

                {
                    OPstream toNeighbour(Pstream::commsTypes::blocking, procPatch.neighbProcNo());

                    toNeighbour << intProcPoints;
                }
            }
        }

        //- receiving
        forAll(mesh.boundaryMesh(), patchI)
        {
            const polyPatch& patch = mesh.boundaryMesh()[patchI];

            if (isA<processorPolyPatch>(patch))
            {
                List<bool> intProcPoints;

                const processorPolyPatch& procPatch =
                        refCast<const processorPolyPatch>(patch);

                {
                    IPstream fromNeighbour(Pstream::commsTypes::blocking, procPatch.neighbProcNo());

                    fromNeighbour >> intProcPoints;
                }

                const labelList& globalPatchPoints = correspPatchPoints_[patchI];

                forAll(intProcPoints, p)
                {
                    if(intProcPoints[p])
                    {
                        const label& patchPointI = globalPatchPoints[p];

                        pointData[patchPointI].pointOnInterface();
                    }
                }
            }
        }
    }
}

// set the region Name
void region::setRegionName
(
    const word& name
)
{
    regionName_ = name;
}

// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
