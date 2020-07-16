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

\*---------------------------------------------------------------------------*/

#include "dsmcParcel.H"
#include "dsmcCloud.H"
#include "meshTools.H"


// * * * * * * * * * * * * * * * Static Data * * * * * * * * * * * * * * * * //

Foam::label Foam::dsmcParcel::TrackedParcel::nDELETED = 0;


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::dsmcParcel::move
(
    dsmcParcel::trackingData& td,
    const scalar trackTime
)
{
    td.switchProcessor = false;
    td.keepParticle = true;

    if (isFree())
    {
        const polyMesh& mesh = td.cloud().pMesh();
        const polyBoundaryMesh& pbMesh = mesh.boundaryMesh();

        if (newParcel() != -1)
        {
            // note: this justifies that freshly inserted parcels should be
            // tracked as if they had passed the boundary face on which they
            // have been inserted in the time step in which they are inserted.
            stepFraction() = td.cloud().rndGen().sample01<scalar>();
            newParcel() = -1;
        }

        //scalar tEnd = (1.0 - stepFraction())*trackTime; // OLD FORMULATION
        label orgCell = cell(); // NEW VINCENT
        scalar tEnd = (1.0 - stepFraction())*td.cloud().deltaTValue(orgCell); // NEW VINCENT
        //const scalar dtMax = tEnd; // OLD FORMULATION

        // For reduced-D cases, the velocity used to track needs to be
        // constrained, but the actual U_ of the parcel must not be
        // altered or used, as it is altered by patch interactions one
        // needs to retain its 3D value for collision purposes.
        vector Utracking = U_;

        while (td.keepParticle && !td.switchProcessor && tEnd > ROOTVSMALL)
        {
            Utracking = U_;

            if (td.cloud().coordSystem().type() == "dsmcCartesian")
            {
                // Apply correction to position for reduced-D cases,
                // but not for axisymmetric cases
                meshTools::constrainToMeshCentre(mesh, position());

                // Apply correction to velocity to constrain tracking for
                // reduced-D cases,  but not for axisymmetric cases
                meshTools::constrainDirection(mesh, mesh.solutionD(), Utracking);
            }

            //- Set the Lagrangian time-step
            //scalar dt = min(dtMax, tEnd); // OLD FORMULATION
            scalar dt = tEnd; // NEW VINCENT

            orgCell = cell(); // NEW VINCENT
            dt *= trackToFace(position() + dt*Utracking, td, true);
            const label destCell = cell(); // NEW VINCENT

            tEnd -= dt;

            stepFraction() = 1.0 - tEnd/td.cloud().deltaTValue(orgCell); // NEW VINCENT
            //stepFraction() = 1.0 - tEnd/trackTime; // OLD FORMULATION

            /*if (destCell != orgCell)
            {
                tEnd *= td.cloud().deltaTValue(destCell)
                    /td.cloud().deltaTValue(orgCell);
            } // NEW VINCENT*/

            //- face tracking info
            if (face() != -1)
            {
                //- measure flux properties
                td.cloud().tracker().trackParcelFaceTransition(*this);
            }

            if (onBoundary() && td.keepParticle)
            {
                if (isA<processorPolyPatch>(pbMesh[patch(face())]))
                {
                    td.switchProcessor = true;
                }

                forAll(td.cloud().boundaries().cyclicBoundaryModels(), c)
                {
                    const labelList& faces = td.cloud().boundaries().cyclicBoundaryModels()[c]->allFaces();

                    if (findIndex(faces, this->face()) != -1)
                    {
                        td.cloud().boundaries().cyclicBoundaryModels()[c]->controlMol(*this, td);
                    }
                }
            }
        }
    }
    else
    {
        //- The stuck particle is considered for desorption
        //  NOTE: there should be a better way to locate the patch than looping through them all
        forAll(td.cloud().boundaries().patchBoundaryModels(), c)
        {
            if (td.cloud().boundaries().patchBoundaryModels()[c]->patchId() == stuck().wallTemperature()[1])
            {
                // then this patch is the "dsmc*Sticking*WallPatch" the particle is stuck on
                td.cloud().boundaries().patchBoundaryModels()[c]->controlParticle(*this, td);
                break;
            }
        }
    }

    return td.keepParticle;
}


bool Foam::dsmcParcel::hitPatch
(
    const polyPatch&,
    trackingData& td,
    const label,
    const scalar,
    const tetIndices&
)
{
    return false;
}


void Foam::dsmcParcel::hitProcessorPatch
(
    const processorPolyPatch&,
    trackingData& td
)
{
    td.switchProcessor = true;
}


void Foam::dsmcParcel::hitWallPatch
(
    const wallPolyPatch& wpp,
    trackingData& td,
    const tetIndices& tetIs
)
{
    //- find which patch has been hit
    label patchIndex = wpp.index();

    const label& patchModelId = td.cloud().boundaries().patchToModelIds()[patchIndex];

    //- apply a boundary model when a molecule collides with this poly patch
    td.cloud().boundaries().patchBoundaryModels()[patchModelId]->controlParticle(*this, td);
}


void Foam::dsmcParcel::hitPatch
(
    const polyPatch& pp,
    trackingData& td
)
{
    //- find which patch has been hit
    label patchIndex = pp.index();

    const label& patchModelId = td.cloud().boundaries().patchToModelIds()[patchIndex];

    //- apply a boundary model when a molecule collides with this poly patch
    td.cloud().boundaries().patchBoundaryModels()[patchModelId]->controlParticle(*this, td);
}


void Foam::dsmcParcel::transformProperties
(
    const tensor& T
)
{
   particle::transformProperties(T);
   U_ = transform(T, U_);
}


void Foam::dsmcParcel::transformProperties
(
    const vector& separation
)
{
    particle::transformProperties(separation);
}

bool Foam::dsmcParcel::relocateStuckParcel
(
    const polyMesh& mesh
)
{
    const polyBoundaryMesh& bMesh = mesh.boundaryMesh();

    // find the closest patch and patch face indices, as this is the patch/face
    // on which this parcel is stuck.
    scalar closestFaceDistance = GREAT;
    label closestPatchi = -1;
    label closestPatchFacei = -1;

    forAll(mesh.cells()[cell()], i)
    {
        const label facei = mesh.cells()[cell()][i];

        // find corresponding boundary patch
        const label patchi = bMesh.whichPatch(facei);

        if (patchi != -1)
        {
            // this is a boundary patch, i.e. a potential candidate for the
            // patch to which parcel is stuck.
            const polyPatch& wpp = bMesh[patchi];

            // patch face index:
            const label patchFacei = wpp.whichFace(facei);

            // calculate distance between parcel position and this patch face:
            pointHit pHit
            (
                wpp[patchFacei].nearestPoint(position(), wpp.points())
            );

            if (pHit.hit())
            {
                const scalar distance = pHit.distance();
                if (distance < closestFaceDistance)
                {
                    // found new closest patch face
                    closestFaceDistance = distance;
                    closestPatchi = patchi;
                    closestPatchFacei = patchFacei;
                }
            }
        }
    }

    // did we find a closest patch face?
    if (closestPatchi != -1)
    {
        // yes, reset stuck parcel information to this patch face
        stuck().wallTemperature()[1] = closestPatchi;
        stuck().wallTemperature()[2] = closestPatchFacei;
        return true;
    }
    // no, indicates fatal error
    return false;
}


// * * * * * * * * * * * * * * * *  IOStream operators * * * * * * * * * * * //

#include "dsmcParcelIO.C"


// ************************************************************************* //
