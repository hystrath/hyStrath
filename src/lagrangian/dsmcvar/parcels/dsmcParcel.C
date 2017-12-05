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

#include <typeinfo>

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::dsmcParcel::move
(
    dsmcParcel::trackingData& td,
    const scalar trackTime
)
{
    td.switchProcessor = false;
    td.keepParticle = true;
    
    const polyMesh& mesh = td.cloud().pMesh();
    const polyBoundaryMesh& pbMesh = mesh.boundaryMesh();
//     Random& rndGen(td.cloud().rndGen());

    if(newParcel() == 1)
    {
        Random& rndGen(td.cloud().rndGen());
        stepFraction() = rndGen.scalar01(); 
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

        // Set the Lagrangian time-step
        scalar dt = min(dtMax, tEnd);
        
        dt *= trackToFace(position() + dt*Utracking, td, true/*, td.cloud().faceTree()*/);

        tEnd -= dt;

        stepFraction() = 1.0 - tEnd/trackTime;

        // - face tracking info
        if( face() != -1 )    //*******
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

    return td.keepParticle;
}


bool Foam::dsmcParcel::move
(
    dsmcParcel::trackingData& td,
    const scalarField& trackTime
)
{
    td.switchProcessor = false;
    td.keepParticle = true;

    const polyMesh& mesh = td.cloud().pMesh();
    const polyBoundaryMesh& pbMesh = mesh.boundaryMesh();
//     Random& rndGen(td.cloud().rndGen());

    if(newParcel() == 1)
    {
        Random& rndGen(td.cloud().rndGen());
        stepFraction() = rndGen.scalar01(); 
        newParcel() = 0;
    }
    
    scalar tEnd = (1.0 - stepFraction())*trackTime[cell()];
    //const scalar dtMax = tEnd; // TODO DELETED VINCENT
                
    // For reduced-D cases, the velocity used to track needs to be
    // constrained, but the actual U_ of the parcel must not be
    // altered or used, as it is altered by patch interactions an
    // needs to retain its 3D value for collision purposes.
    vector Utracking = U_;
    //label  i=0; 

    while (td.keepParticle && !td.switchProcessor && tEnd > ROOTVSMALL)
    {
        //i++;
        const scalar originalDeltaT = trackTime[cell()]; // NEW VINCENT
        
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

        // Set the Lagrangian time-step
        scalar dt = tEnd; //min(dtMax, tEnd); // TODO dtMax not adjusted when changing cell, should it need to be
        
        //if (i==1) Info << "\n" << endl;
        //Info << "B"<< cell() <<": tEnd" << tab << tEnd << tab << "stepFraction()" << tab << stepFraction() << endl;

        dt *= trackToFace(position() + dt*Utracking, td, true/*, td.cloud().faceTree()*/);
        
        //scalar timeFractionPerformed = trackToFace(position() + dt*Utracking, td, true); // NEW VINCENT
        //dt *= timeFractionPerformed;
        //Info << "I: stepFraction()" << tab << timeFractionPerformed << endl;

        tEnd -= dt;

        stepFraction() = 1.0 - tEnd/originalDeltaT;
        
        //Info << "I: tEnd" << tab << tEnd << tab << "stepFraction()" << tab << stepFraction() << endl;
        
        if(tEnd > ROOTVSMALL) // NEW VINCENT
        {
            const scalar destinationDeltaT = trackTime[cell()];
            
            tEnd *= destinationDeltaT/originalDeltaT; 
            
            //stepFraction() = 1.0 - tEnd/destinationDeltaT; // NOTE VINCENT: USELESS CAUSE EQUIVALENT TO stepFraction() = 1.0 - tEnd,NOT_RESCALED/originalDeltaT;
        }
        //Info << "A"<<cell()<<": tEnd" << tab << tEnd << tab << "stepFraction()" << tab << stepFraction() << endl;
        //Info << "(" << td.keepParticle << !td.switchProcessor << ")" << endl;
               
        // - face tracking info
        if( face() != -1 )    //*******
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
    //-find which patch has been hit
    label patchIndex = wpp.index();

    const label& patchModelId = td.cloud().boundaries().patchToModelIds()[patchIndex];

    // apply a boundary model when a molecule collides with this poly patch
    td.cloud().boundaries().patchBoundaryModels()[patchModelId]->controlParticle(*this, td);
}

void Foam::dsmcParcel::hitPatch
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
void Foam::dsmcParcel::transformProperties
(
    const tensor& T
)
{
   particle::transformProperties(T);
   U_ = transform(T, U_);
}


 //template<class ParcelType>
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
