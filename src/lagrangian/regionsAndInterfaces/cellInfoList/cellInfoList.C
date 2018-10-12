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
    cellInfoList

Description

\*----------------------------------------------------------------------------*/

#include "cellInfoList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Null constructor
cellInfoList::cellInfoList()
:
    List<cellInform> ()
{}

// constructor from list
cellInfoList::cellInfoList
(
    const List<cellInform>& cellInfos 
)
:
    List<cellInform> (cellInfos)
{}

// constructor from components
cellInfoList::cellInfoList
(
    const polyMesh& mesh,
    const faceInfoList& faceData
)
:
    List<cellInform> (mesh.nCells())
{
    setCells();

    setBoundaryCells(mesh, faceData);
}



// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

// Ostream& operator<<(Ostream& os, const cellInfoList& p)
// {
//     if (os.format() == IOstream::ASCII)
//     {
//         return os
//             << p.currentPointLabel() << p.originPoint() <<  p.Rsqr();
//     }
//     else
//     {
//         return os
//             <<  p.currentPointLabel() << p.originPoint() <<  p.Rsqr();
//     }
// }
// 
// Istream& operator>>(Istream& is, cellInfoList& p)
// {
//     return is >> p.currentPointLabel_ >> p.originPoint_ >>  p.Rsqr_;
// }


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

cellInfoList::~cellInfoList()
{}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void cellInfoList::setBoundaryCells
(
    const polyMesh& mesh,
    const faceInfoList& faceData
)
{
    forAll(*this, c)
    {
        const labelList& facesI = mesh.cells()[c];

        operator[](c).checkCellBoundary(faceData, facesI);
    }
}

void cellInfoList::setCells()
{
    forAll(*this, c)
    {
//         const labelList& facesI = mesh.cells()[c];
        operator[](c).setCellInfo(label(c));
    }
}

void cellInfoList::cellsInZone(const labelList& cellZone)
{
    forAll(cellZone, c)
    {
        operator[](cellZone[c]).cellInZone();
    }
}


// void cellInfoList::cellsInPrefZone(const labelList& cellZone)
// {
//     forAll(cellZone, c)
//     {
//         operator[](cellZone[c]).cellInPrefZone();
//     }
// }

void cellInfoList::resetAcceptedCells(const labelList& cellZone)
{
    forAll(cellZone, c)
    {
        operator[](cellZone[c]).resetAcceptedCell();
    }
}

void cellInfoList::cellsOutOfZone(const labelList& cellZone)
{
    forAll(cellZone, c)
    {
        operator[](cellZone[c]).cellOutOfZone();
    }
}

// void cellInfoList::cellsOutOfPrefZone(const labelList& cellZone)
// {
//     forAll(cellZone, c)
//     {
//         operator[](cellZone[c]).cellOutOfPrefZone();
//     }
// }


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //





// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
