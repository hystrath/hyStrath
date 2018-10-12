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

#include "pdParcel.H"
#include "pdCloud.H"
#include "meshTools.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::pdParcel::move
(
    pdParcel::trackingData& td,
    const scalar trackTime
)
{
    td.switchProcessor = false;
    td.keepParticle = true;

    const polyMesh& mesh = td.cloud().pMesh();
    const polyBoundaryMesh& pbMesh = mesh.boundaryMesh();

    //if the parcel is new
    if(newParcel() == 1)
    {
        Random& rndGen(td.cloud().rndGen());
        stepFraction() = td.cloud().rndGen().sample01<scalar>(); 

        //roll back velocity
        U_ -= 0.5*A_*trackTime;

        newParcel() = 0;
    }

    // First leapfrog velocity adjust part, required before tracking+force part
    if(td.part() == 0)
    {
        U_ += A_*trackTime;
        //reduce(U_,sumOp<vector>());
    }
    // Leapfrog tracking part
    else if (td.part() == 1)
    {
        scalar tEnd = (1.0 - stepFraction())*trackTime;
        scalar dtMax = tEnd;

        // For reduced-D cases, the velocity used to track needs to be
        // constrained, but the actual U_ of the parcel must not be
        // altered or used, as it is altered by patch interactions an
        // needs to retain its 3D value for collision purposes.
        vector Utracking = U_;

        //moves the particle until it has either moved to x_k+1, moved processor or died.
        while (td.keepParticle && !td.switchProcessor && tEnd > ROOTVSMALL)
        {
            // Apply correction to position for reduced-D cases
            meshTools::constrainToMeshCentre(mesh, position());

            Utracking = U_;

            // Apply correction to velocity to constrain tracking for
            // reduced-D cases

            meshTools::constrainDirection(mesh, mesh.solutionD(), Utracking);

            // Set the Lagrangian time-step
            scalar dt = min(dtMax, tEnd);

            //track particle to a face
            dt *= trackToFace(position() + dt*Utracking, td, true);

            //remove time taken to get to a face
            tEnd -= dt;

            //Info << "   tEnd " << tEnd << endl;
            stepFraction() = 1.0 - tEnd/trackTime;

            // - face tracking info
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
    else if (td.part() == 2)
    {
        //Roll velocity backward by half a step
        U_ -= 0.5*A_*trackTime;
        //reduce(U_,sumOp<vector>());
    }
    else if (td.part() == 3)
    {
        //Roll velocity forward by half a step
        U_ += 0.5*A_*trackTime;
        //reduce(U_,sumOp<vector>());
    }
    //Error catch
    else
    {
        FatalErrorIn("pdParcel::move(trackingData&, const scalar)") << nl
            << td.part() << " is an invalid part of the integration method."
            << abort(FatalError);
    }
    return td.keepParticle;
}

bool Foam::pdParcel::hitPatch
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

void Foam::pdParcel::hitProcessorPatch
(
    const processorPolyPatch&,
    trackingData& td
)
{
    td.switchProcessor = true;
}

void Foam::pdParcel::hitWallPatch
(
    const wallPolyPatch& wpp,
    trackingData& td,
    const tetIndices& tetIs
)
{
    //-find which patch has been hit
    label patchIndex = wpp.index();

    const label& patchModelId = td.cloud().boundaries().patchToModelIds()[patchIndex];

    // apply a boundary model when a molecule collides with this poly patch
    td.cloud().boundaries().patchBoundaryModels()[patchModelId]->controlParticle(*this, td);
}

void Foam::pdParcel::hitPatch
(
    const polyPatch& pp,
    trackingData& td
)
{
    //-find which patch has been hit
    label patchIndex = pp.index();

    const label& patchModelId = td.cloud().boundaries().patchToModelIds()[patchIndex];

    // apply a boundary model when a molecule collides with this poly patch
    td.cloud().boundaries().patchBoundaryModels()[patchModelId]->controlParticle(*this, td);
}


 //template<class ParcelType>
void Foam::pdParcel::transformProperties
(
    const tensor& T
)
{
   particle::transformProperties(T);
   U_ = transform(T, U_);
}


 //template<class ParcelType>
void Foam::pdParcel::transformProperties
(
    const vector& separation
)
{
  particle::transformProperties(separation);
}


/*Foam::pdParcel::~pdParcel()
{
    //Info << "pdParcel Destructor" << endl;
}*/


// * * * * * * * * * * * * * * * *  IOStream operators * * * * * * * * * * * //

#include "pdParcelIO.C"

// ************************************************************************* //
