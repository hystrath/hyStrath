/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2009-2010 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

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

Foam::scalar Foam::dsmcParcel::TrackedParcel::D_EFF = 0.0;


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::dsmcParcel::move
(
    dsmcParcel::trackingData& td,
    const scalar trackTime
)
{
    td.switchProcessor = false;
    td.keepParticle = true;

    if(isFree())
    {
        const polyMesh& mesh = td.cloud().pMesh();
        const polyBoundaryMesh& pbMesh = mesh.boundaryMesh();

        if(newParcel() == 1)
        {
            stepFraction() = td.cloud().rndGen().scalar01(); 
            newParcel() = 0;
        }
        
        
        scalar tEnd = (1.0 - stepFraction())*trackTime;
        const scalar dtMax = tEnd;
                    
        // For reduced-D cases, the velocity used to track needs to be
        // constrained, but the actual U_ of the parcel must not be
        // altered or used, as it is altered by patch interactions an
        // needs to retain its 3D value for collision purposes.
        vector Utracking = U_;
            
        while (td.keepParticle && !td.switchProcessor && tEnd > ROOTVSMALL)
        {
            Utracking = U_;

            if(!td.cloud().axisymmetric())
            {
                // Apply correction to position for reduced-D cases, 
                // but not axisymmetric cases
                meshTools::constrainToMeshCentre(mesh, position());
                
                // Apply correction to velocity to constrain tracking for
                // reduced-D cases,  but not axisymmetric cases
                meshTools::constrainDirection(mesh, mesh.solutionD(), Utracking);
            }

            //- Set the Lagrangian time-step
            scalar dt = min(dtMax, tEnd);

            dt *= trackToFace(position() + dt*Utracking, td, true/*, td.cloud().faceTree()*/);

            tEnd -= dt;

            stepFraction() = 1.0 - tEnd/trackTime;
                
            //- face tracking info
            if( face() != -1 )
            {
                //--  measure flux properties
                td.cloud().tracker().updateFields
                (
                    *this
                );
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

                    if(findIndex(faces, this->face()) != -1)
                    {
                        td.cloud().boundaries().cyclicBoundaryModels()[c]->controlMol(*this, td);
                    }
                }
            }
        }
    }
    else // NEW VINCENT
    {
        //- The stuck particle is considered for desorption
        // NOTE: there should be a better way to locate the patch than looping through them all
        forAll(td.cloud().boundaries().patchBoundaryModels(), c)
        {
            if(td.cloud().boundaries().patchBoundaryModels()[c]->patchId() == stuck().wallTemperature()[1])
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


// * * * * * * * * * * * * * * * *  IOStream operators * * * * * * * * * * * //

#include "dsmcParcelIO.C"

// ************************************************************************* //
