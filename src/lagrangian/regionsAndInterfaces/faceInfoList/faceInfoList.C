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
    faceInfoList

Description

\*----------------------------------------------------------------------------*/

#include "faceInfoList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Null constructor
faceInfoList::faceInfoList()
:
    List<faceInfo> ()
{}

// constructor from list
faceInfoList::faceInfoList
(
    const List<faceInfo>& faceInfos 
)
:
    List<faceInfo> (faceInfos)
{}

// constructor from components
faceInfoList::faceInfoList
(
    const polyMesh& mesh
)
:
    List<faceInfo> (mesh.nFaces())
{
    setFaces();

    setBoundaryFaces(mesh);
}



// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

// Ostream& operator<<(Ostream& os, const faceInfoList& p)
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
// Istream& operator>>(Istream& is, faceInfoList& p)
// {
//     return is >> p.currentPointLabel_ >> p.originPoint_ >>  p.Rsqr_;
// }


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

faceInfoList::~faceInfoList()
{}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void faceInfoList::setBoundaryFaces( const polyMesh& mesh )
{
    forAll( mesh.boundaryMesh(), patchI )
    {
        const polyPatch& patch = mesh.boundaryMesh()[patchI];

        label startFaceI = 0;

//         Info << " Patch start: " << patch.start() << endl;

        for(label i = 0; i < patch.size(); i++)
        {
            label patchFaceI = i + startFaceI;

            label meshFaceI = patch.start() + patchFaceI;

//             Info << " Face label: " << meshFaceI << endl;

            // determining the boundary faces

            operator[](meshFaceI).faceOnBoundary();

            if (isA<processorPolyPatch>(patch))
            {
                operator[](meshFaceI).faceOnProcBoundary();
            }

            if (isA<cyclicPolyPatch>(patch))
            {
                operator[](meshFaceI).faceOnCyclicBoundary();
            }
        }
    }
}


void faceInfoList::setFaces()
{
    forAll(*this, f)
    {
//         const labelList& pointsI = mesh.faces()[f];

        operator[](f).setFaceInfo(label(f));
    }
}


void faceInfoList::facesOnInterface( const labelList& faces )
{
    forAll(faces, f)
    {
        operator[](faces[f]).faceOnInterface();
    }
}


void faceInfoList::removeAcceptedFaces()
{
    forAll(*this, f)
    {
        operator[](f).resetAcceptedFace();
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //





// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
