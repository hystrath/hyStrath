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
    dsmcFaceTracker

Description

\*----------------------------------------------------------------------------*/

#include "dsmcFaceTracker.H"
#include "processorPolyPatch.H"
#include "cyclicPolyPatch.H"
#include "wallPolyPatch.H"
#include "dsmcCloud.H"

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


// Construct from mesh and cloud 
dsmcFaceTracker::dsmcFaceTracker
(
    const polyMesh& mesh,
    dsmcCloud& cloud
)
:
    mesh_(refCast<const fvMesh>(mesh)),
    cloud_(cloud)
{}


//- Construct from mesh, cloud and boolean (dsmcFoam)
dsmcFaceTracker::dsmcFaceTracker
(
    const polyMesh& mesh,
    dsmcCloud& cloud,
    const bool& dummy
)
:
    mesh_(refCast<const fvMesh>(mesh)),
    cloud_(cloud),
    parcelIdFlux_(cloud_.typeIdList().size()),
    massIdFlux_(cloud_.typeIdList().size())

{
    forAll(parcelIdFlux_, i)
    {
        parcelIdFlux_[i].setSize(mesh_.nFaces(), 0.0);
        massIdFlux_[i].setSize(mesh_.nFaces(), 0.0);
    }
}
// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

dsmcFaceTracker::~dsmcFaceTracker()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void dsmcFaceTracker::clean()
{
    //- clean geometric fields

    forAll(parcelIdFlux_, i)
    {
        parcelIdFlux_[i] = scalar(0.0);
        massIdFlux_[i] = scalar(0.0);
    }
}


//reset size of fields after mesh has changed
void dsmcFaceTracker::reset()
{
    forAll(parcelIdFlux_, i)
    {
        parcelIdFlux_[i].setSize(mesh_.nFaces(), 0.0);
        massIdFlux_[i].setSize(mesh_.nFaces(), 0.0);
    }   
}


void dsmcFaceTracker::updateFields
(
    dsmcParcel& p
)
{
    const label& crossedFace = p.face();
    const label& typeId = p.typeId();
    const dsmcParcel::constantProperties& constProp = cloud_.constProps(typeId);
    const scalar& mass = constProp.mass();
    const vector& U = p.U();
//     const vector mom = p.U()*mass;

    //- check which patch was hit
    const label& patchId = mesh_.boundaryMesh().whichPatch(crossedFace);

    //- direction of dsmcParcel trajectory with respect to the face normal
    scalar sgn = sign( U & mesh_.faceAreas()[crossedFace] ) * 1.0;

    //- geometric fields

    if(patchId != -1) //- boundary face
    {
        const polyPatch& patch = mesh_.boundaryMesh()[patchId];

        const label faceIndex = crossedFace - patch.start();

        //- correct cyclic patches
        if (isA<cyclicPolyPatch>(patch))
        {
            label coupledFace = refCast<const cyclicPolyPatch>
            (
                patch
            ).neighbPatch().start() + faceIndex;

            parcelIdFlux_[typeId][coupledFace] += 1.0;
            massIdFlux_[typeId][coupledFace] += mass;
        }


        // molecular properties are appended to the face of the leaving processor only.
        // normal vector points out from the domain.
        if (isA<processorPolyPatch>(patch))
        {
            parcelIdFlux_[typeId][crossedFace] += sgn*1.0;
            massIdFlux_[typeId][crossedFace] += sgn*mass;
        }
    }
    else //- internal face
    {
        //- properties
        parcelIdFlux_[typeId][crossedFace] += sgn*1.0;
        massIdFlux_[typeId][crossedFace] += sgn*mass;
    }
}

// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

}  // End namespace Foam

// ************************************************************************* //
