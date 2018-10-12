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
    cellRegions

Description

\*----------------------------------------------------------------------------*/

#include "cellRegions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //





// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

cellRegions::~cellRegions()
{}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

//- construct null
// cellRegions::cellRegions()
// :
//     List<subLayer> ()
// {}

//- construct from components
// cellRegions::cellRegions
// (
//     const List<subLayer>& subLayerList,
// )
// :
//     List<subLayer> (subLayerList),
// {}


//- Construct from components: polyMesh, list of pointInfo,
//   list of faces (interface), list of cells (zone) and
//   a list of sub layer thicknesses.

cellRegions::cellRegions
(
    const polyMesh& mesh,
    pointInfoList& pointData,
    faceInfoList& faceData,
    cellInfoList& cellData,
    const labelList& zone,
    const scalar& radius
)
:
    List<region>(),
    Rsqr_(sqr(radius)),
    cellProcAddressing_()
{
    label nCells = mesh.nCells();

    if (Pstream::parRun())
    {
        reduce(nCells, sumOp<label>());
    }

    if (Pstream::parRun())
    {
        labelIOList cellProcAddressing
        (
            IOobject
            (
                "cellProcAddressing",
                mesh.facesInstance(),
                mesh.meshSubDir,
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        );

        cellProcAddressing_ = cellProcAddressing;

        Info << "cellProcAddressing: " << cellProcAddressing_ << endl;
    }

    List<region> regions(nCells);

    forAll(regions, r)
    {
        Info << "Building region: " << r << endl;

        DynamicList<label> interface(0);

        if (Pstream::parRun())
        {
            const label& localCell = findIndex(cellProcAddressing_, r);

            if(localCell != -1)
            {
                const labelList& faceList = mesh.cells()[localCell];

                forAll(faceList, f)
                {
                    interface.append(faceList[f]);
                }
            }
        }
        else
        {
            const labelList& faceList = mesh.cells()[r];

            forAll(faceList, f)
            {
                interface.append(faceList[f]);
            }
        }

        std::ostringstream oss;
        oss << r;
        std::string s(oss.str());

        const word regionName = "region_" + word(s);

        regions[r].setRegion
        (
            mesh,
            zone,
            interface.shrink(),
            Rsqr_,
            regionName,
            pointData,
            faceData,
            cellData
        );
    }

    transfer(regions);
}


// cells belonging to all the regions
// void cellRegions::setRegionsCells
// (
//     const labelList& cellsInRegion
// )
// {
//     DynamicList<label> regionsCells(regionsCells_);
// 
//     forAll(cellsInRegion, c)
//     {
//         regionsCells.append(cellsInRegion[c]);
//     }
// 
//     regionsCells_ = regionsCells.shrink();
// }

// cells remaining within the zone 
// i.e. cells that are in the zone but not part of the regions
// void cellRegions::setRemainingCells
// (
//     const labelList& regionCells
// )
// {
//     DynamicList<label> newZone(0);
// 
//     forAll(remainingCells_, c)
//     {
//         if(findIndex(regionCells, remainingCells_[c]) == -1)
//         {
//             newZone.append(remainingCells_[c]);
//         }
//     }
// 
//     remainingCells_ = newZone.shrink();
// }







// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
