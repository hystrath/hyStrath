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
    interface

Description

\*----------------------------------------------------------------------------*/

#include "interface.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

interface::interface
(
    const polyMesh& mesh,
    faceInfoList& faceData,
    cellInfoList& cellData,
    const labelList& cellZoneA,
    const word& zoneNameA,
    const labelList& cellZoneB,
    const word& zoneNameB
)
:
    mesh_(mesh),
    faceData_(faceData),
    cellData_(cellData),
    cellZoneA_(cellZoneA),
    zoneNameA_(zoneNameA),
    cellZoneB_(cellZoneB),
    zoneNameB_(zoneNameB),
    boundaryCells_(0),
    interface_(0)
{
    setInterface();
    setName();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

interface::~interface()
{}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void interface::setInterface()
{
    // an interface cannot be built if between two adjacent zones, one zone has no cells.
    if(
        (cellZoneA_.size() > 0) &&
        (cellZoneB_.size() > 0)
       )
    {
        setInternalInterface();
    }

    boundaryCells_.shrink();

    handleBoundaries();

//     syncInterfaceFaces();

    interface_.shrink();
//     faceData_.facesOnInterface(interface_); // be careful of this -- i don't agree

    // reset data
    cellData_.resetAcceptedCells(cellZoneA_);
    faceData_.removeAcceptedFaces();
}

void interface::setInternalInterface()
{
    // loop through all cells within cellZoneA and determine the interface with cellZoneB
    forAll(cellZoneA_, cellA)
    {
        const label& cellI = cellZoneA_[cellA];

        //accept the cell as being tested
        cellData_[cellI].acceptCell();

        const labelList& neighbourCells = mesh_.cellCells()[cellI];

        forAll(neighbourCells, cell)
        {
            // neighbour cell
            const label& cellN = neighbourCells[cell];

            if(!cellData_[cellN].isAcceptedCell())
            {
                //if neighbour cell is in cellZoneB
                if(findIndex(cellZoneB_, cellN) != -1)
                {
                    // append the common face to the interface
                    findInterfaceFace(cellI, cellN);
                }
            }
        }

        //handle cells on the boundary (if parallel processing)
        if(Pstream::parRun())
        {
            if(cellData_[cellI].isBoundaryCell())
            {
                boundaryCells_.append(cellI);
            }
        }
    }
}

inline void interface::findInterfaceFace
(
    const label& cellI,
    const label& cellN
)
{
    const labelList& facesCellI = mesh_.cells()[cellI];
    const labelList& facesCellN = mesh_.cells()[cellN];

    forAll(facesCellI, faceI)
    {
        forAll(facesCellN, faceN)
        {
            if(facesCellI[faceI] == facesCellN[faceN])
            {
                interface_.append(facesCellI[faceI]);
            }
        }
    }
}


void interface::handleBoundaries()
{
    if(Pstream::parRun())
    {
        // for all cells that have been visited at the boundary of the processor domain 
        // accept the faces at the processor boundary.
        forAll(boundaryCells_, cell)
        {
            const label& cellI = boundaryCells_[cell];

            const labelList& connectedFaces = mesh_.cells()[cellI];

            // loop through all faces
            forAll(connectedFaces, face)
            {
                const label& faceI = connectedFaces[face];

                // if face is on a processor boundary
                if(faceData_[faceI].isProcBoundaryFace())
                {
                    faceData_[faceI].acceptFace();
                }
            }
        }


        //sending data to neighbouring processors
        forAll(mesh_.boundaryMesh(), patchI)
        {
            const polyPatch& patch = mesh_.boundaryMesh()[patchI];

            if (isA<processorPolyPatch>(patch))
            {
                List<bool> patchFaces(patch.size());

                forAll(patchFaces, i)
                {
                    label meshFaceI = patch.start() + label(i);

                    patchFaces[i] = faceData_[meshFaceI].isAcceptedFace();
                }

                // sending information across neighbouring processors.
                const processorPolyPatch& procPatch =
                        refCast<const processorPolyPatch>(patch);

                {
                    OPstream toNeighbour(Pstream::commsTypes::blocking, procPatch.neighbProcNo());

                    toNeighbour << patchFaces;
                }
            }
        }

        faceData_.removeAcceptedFaces();

        //receiving data from neighbouring processors
        forAll(mesh_.boundaryMesh(), patchI)
        {
            const polyPatch& patch = mesh_.boundaryMesh()[patchI];

            if (isA<processorPolyPatch>(patch))
            {
                const processorPolyPatch& procPatch =
                    refCast<const processorPolyPatch>(patch);

                List<bool> patchFaces;

                {
                    IPstream fromNeighbour(Pstream::commsTypes::blocking, procPatch.neighbProcNo());

                    fromNeighbour >> patchFaces;
                }

                forAll(patchFaces, i)
                {
                    // if face is accepted check cell within zone
                    if(patchFaces[i])
                    {
                        label meshFaceI = patch.start() + i;

                        const label& cellP = mesh_.faceOwner()[meshFaceI];

                        if(!cellData_[cellP].isAcceptedCell())
                        {
                            //if cell is in cellZoneB
                            if(findIndex(cellZoneB_, cellP) != -1)
                            {
                                interface_.append(meshFaceI);

                                faceData_[meshFaceI].isAcceptedFace();
                            }
                        }
                    }
                }
            }
        }
    }
}


void interface::syncInterfaceFaces()
{
    if(Pstream::parRun())
    {
        //- sending
        forAll(mesh_.boundaryMesh(), patchI)
        {
            const polyPatch& patch = mesh_.boundaryMesh()[patchI];

            if (isA<processorPolyPatch>(patch))
            {
                List<bool> patchFaces(patch.size());

                forAll(patchFaces, i)
                {
                    label meshFaceI = patch.start() + label(i);

                    patchFaces[i] = faceData_[meshFaceI].isAcceptedFace();
                }

                // sending information across neighbouring processors.
                const processorPolyPatch& procPatch =
                        refCast<const processorPolyPatch>(patch);

                {
                    OPstream toNeighbour(Pstream::commsTypes::blocking, procPatch.neighbProcNo());

                    toNeighbour << patchFaces;
                }
            }
        }

        //receiving data from neighbouring processors
        forAll(mesh_.boundaryMesh(), patchI)
        {
            const polyPatch& patch = mesh_.boundaryMesh()[patchI];

            if (isA<processorPolyPatch>(patch))
            {
                const processorPolyPatch& procPatch =
                    refCast<const processorPolyPatch>(patch);

                List<bool> patchFaces;

                {
                    IPstream fromNeighbour(Pstream::commsTypes::blocking, procPatch.neighbProcNo());

                    fromNeighbour >> patchFaces;
                }

                forAll(patchFaces, i)
                {
                    // if face is accepted check cell within zone
                    if(patchFaces[i])
                    {
                        label meshFaceI = patch.start() + i;

                        interface_.append(meshFaceI);
                    }
                }
            }
        }

    }
}

void interface::setName()
{
    interfaceName_ = zoneNameA_ + "_" + zoneNameB_;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const labelList& interface::interfaceFaces() const
{
    return interface_;
}


const word& interface::interfaceName() const
{
    return interfaceName_;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
