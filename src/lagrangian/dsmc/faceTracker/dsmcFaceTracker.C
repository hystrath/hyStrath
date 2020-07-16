/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016-2020 hyStrath
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of hyStrath, a derivative work of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

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


void dsmcFaceTracker::trackParcelFaceTransition
(
    const dsmcParcel& p
)
{
    // parcel properties:
    const label& typeId = p.typeId();
    const vector& U = p.U();
    const scalar& RWF = p.RWF();
    const label& crossedFaceI = p.face();

    trackFaceTransition(typeId, U, RWF, crossedFaceI);
}


void dsmcFaceTracker::trackFaceTransition
(
    const label& typeId,
    const vector& U,
    const scalar& RWF,
    const label& crossedFaceI
)
{
    // Note: We have to use the parcels RWF in the following accumulations, the
    // reason for this being that the step in which parcels might be cloned or
    // deleted (e.g. in axisymmetric simulations) is carried out _after_ the
    // parcel movement step is completed. Therefore we have to use the current
    // parcel RWF here. Obviously this also has to be considered in all further
    // calculations that are based on parcelIdFlux and massIdFlux.
    const dsmcParcel::constantProperties& constProp = cloud_.constProps(typeId);
    const scalar& mass = constProp.mass();
    // TODO: This is the place to add momentum (m*U), kinetic energy, etc. if
    // these shall also be counted.

    //- check which patch was hit
    const label& patchId = mesh_.boundaryMesh().whichPatch(crossedFaceI);

    //- direction of dsmcParcel trajectory with respect to the face normal
    scalar sgn = sign( U & mesh_.faceAreas()[crossedFaceI] ) * 1.0;

    // check patch type and count accordingly
    if(patchId != -1) // face belongs to boundary patch
    {
        const polyPatch& patch = mesh_.boundaryMesh()[patchId];
        const label faceIndex = crossedFaceI - patch.start();

        // check for the most likely case first. In general that should mean
        // processor patches.

        // processor patches:
        //   Molecular properties are appended to the face of the leaving
        //   processor only. Normal vector points out from the domain.
        if (isA<processorPolyPatch>(patch))
        {
            parcelIdFlux_[typeId][crossedFaceI] += sgn*RWF;
            massIdFlux_[typeId][crossedFaceI] += sgn*RWF*mass;
        }

        // cyclic patches:
        else if (isA<cyclicPolyPatch>(patch))
        {
            label coupledFace = refCast<const cyclicPolyPatch>
            (
                patch
            ).neighbPatch().start() + faceIndex;

            parcelIdFlux_[typeId][coupledFace] += RWF;
            massIdFlux_[typeId][coupledFace] += RWF*mass;
        }

        // all other boundary patches:
        else
        {
            parcelIdFlux_[typeId][crossedFaceI] += sgn*RWF;
            massIdFlux_[typeId][crossedFaceI] += sgn*RWF*mass;
        }
    }
    else // internal face
    {
        parcelIdFlux_[typeId][crossedFaceI] += sgn*RWF;
        massIdFlux_[typeId][crossedFaceI] += sgn*RWF*mass;
    }
}

// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

}  // End namespace Foam

// ************************************************************************* //
