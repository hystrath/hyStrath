/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2008-2009 OpenCFD Ltd.
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


\*---------------------------------------------------------------------------*/

#include "polyMoleculeCloud.H"
#include "polyMolecule.H"
#include "Time.H"
#include "cyclicPolyPatch.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::tensor Foam::polyMolecule::rotationTensorX(scalar phi) const
{
    return tensor
    (
        1, 0, 0,
        0, Foam::cos(phi), -Foam::sin(phi),
        0, Foam::sin(phi), Foam::cos(phi)
    );
}


Foam::tensor Foam::polyMolecule::rotationTensorY(scalar phi) const
{
    return tensor
    (
        Foam::cos(phi), 0, Foam::sin(phi),
        0, 1, 0,
        -Foam::sin(phi), 0, Foam::cos(phi)
    );
}


Foam::tensor Foam::polyMolecule::rotationTensorZ(scalar phi) const
{
    return tensor
    (
        Foam::cos(phi), -Foam::sin(phi), 0,
        Foam::sin(phi), Foam::cos(phi), 0,
        0, 0, 1
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


bool Foam::polyMolecule::move
(
    polyMolecule::trackingData& td,
    const scalar& trackTime
)
{
    td.switchProcessor = false;
    td.keepParticle = true;

    if (special_ != SPECIAL_FROZEN)
    {
        scalar tEnd = (1.0 - stepFraction())*trackTime;
        scalar dtMax = tEnd;

        while (td.keepParticle && !td.switchProcessor && tEnd > ROOTVSMALL)
        {
            // set the lagrangian time-step
            scalar dt = min(dtMax, tEnd);

            dt *= trackToFace(position() + dt*v_, td, false);

            //- face tracking info
            if( this->face() != -1 )   
            {
                //--  monitoring flux properties
                td.cloud().tracker().updateFields
                (
                    *this
                );
            }

            tEnd -= dt;
            stepFraction() = 1.0 - tEnd/trackTime;
        }
    }

    return td.keepParticle;
}

void Foam::polyMolecule::setAsReferred()
{
    special_ = 1;
}

void Foam::polyMolecule::updateHalfVelocity
(
//     const constantProperties& constProps,
    const constantMoleculeProperties& cP,       
    const scalar& trackTime
)
{
    v_ += 0.5*trackTime*a_;

    pi_ += 0.5*trackTime*tau_;

    if (cP.pointMolecule(id_))
    {
        tau_ = vector::zero;

        pi_ = vector::zero;
    }

    if (cP.linearMolecule(id_))
    {
        tau_.x() = 0.0;

        pi_.x() = 0.0;
    }
}

void Foam::polyMolecule::updateAcceleration
(
//     const constantProperties& constProps
    const constantMoleculeProperties& cP    
)
{
//     scalar m = constProps.mass();
    scalar m = cP.mass(id_);

    forAll(siteForces_, s)
    {
        const vector& f = siteForces_[s];

        a_ += f/m;

        tau_ += (cP.siteRefPositions()[id_][s] ^ (Q_.T() & f));
    }
}

void Foam::polyMolecule::updateAfterMove
(
//     const constantProperties& constProps,
    const constantMoleculeProperties& cP,    
    const scalar& trackTime
)
{   
    if (!cP.pointMolecule(id_))
    {
        const diagTensor& momentOfInertia(cP.momentOfInertia(id_));

        tensor R;

        if (!cP.linearMolecule(id_))
        {
            R = rotationTensorX(0.5*trackTime*pi_.x()/momentOfInertia.xx());
            pi_ = pi_ & R;
            Q_ = Q_ & R;
        }

        R = rotationTensorY(0.5*trackTime*pi_.y()/momentOfInertia.yy());
        pi_ = pi_ & R;
        Q_ = Q_ & R;

        R = rotationTensorZ(trackTime*pi_.z()/momentOfInertia.zz());
        pi_ = pi_ & R;
        Q_ = Q_ & R;

        R = rotationTensorY(0.5*trackTime*pi_.y()/momentOfInertia.yy());
        pi_ = pi_ & R;
        Q_ = Q_ & R;

        if (!cP.linearMolecule(id_))
        {
            R = rotationTensorX(0.5*trackTime*pi_.x()/momentOfInertia.xx());
            pi_ = pi_ & R;
            Q_ = Q_ & R;
        }
    }

    setSitePositions(cP);
}


void Foam::polyMolecule::transformProperties(const tensor& T)
{
    particle::transformProperties(T);
    
    Q_ = T & Q_;

    v_ = transform(T, v_);

    a_ = transform(T, a_);

    pi_ = Q_.T() & transform(T, Q_ & pi_);

    tau_ = Q_.T() & transform(T, Q_ & tau_);

    rf_ = transform(T, rf_);

    sitePositions_ = position_ + (T & (sitePositions_ - position_));

    siteForces_ = T & siteForces_;    
}


void Foam::polyMolecule::transformProperties(const vector& separation)
{
    particle::transformProperties(separation);

    if (special_ == SPECIAL_TETHERED)
    {
        specialPosition_ += separation;
    }

    sitePositions_ = sitePositions_ + separation;
}

/*
void Foam::polyMolecule::setSitePositions(const constantProperties& constProps)
{
    forAll(constProps.sites(), s)
    {
        sitePositions_[s] = position_ + (Q_ & constProps.sites()[s].siteReferencePosition());
    }
}*/

void Foam::polyMolecule::setSitePositions(const constantMoleculeProperties& cP)
{
    forAll(cP.siteRefPositions()[id_], s)
    {
        sitePositions_[s] = position_ + (Q_ & cP.siteRefPositions()[id_][s]);
    }
}

void Foam::polyMolecule::setSiteSizes(label size)
{
    sitePositions_.setSize(size);
    siteForces_.setSize(size);
}

bool Foam::polyMolecule::hitPatch
(
    const polyPatch&,
    trackingData&,
    const label,
    const scalar,
    const tetIndices&
)
{
    return false;
}

void Foam::polyMolecule::hitProcessorPatch
(
    const processorPolyPatch&,
    trackingData& td
)
{
    td.switchProcessor = true;
}

void Foam::polyMolecule::hitWallPatch
(
    const wallPolyPatch& wpp,
    trackingData& td,
    const tetIndices& tetIs
)
{
    // Use of the normal from tetIs is not required as
    // hasWallImpactDistance for a moleculeCloud is false.
    vector nw = normal();
    nw /= mag(nw);

    scalar vn = v_ & nw;

    // Specular reflection
    if (vn > 0)
    {
        v_ -= 2*vn*nw;
    }
}

void Foam::polyMolecule::hitPatch
(
    const polyPatch& pp,
    trackingData& td
)
{
    //-find which patch has been hit
    label patchIndex = pp.index();

    const label& patchModelId = td.cloud().boundaries().patchToModelIds()[patchIndex];

    // apply a boundary model when a molecule collides with this poly patch
    td.cloud().boundaries().patchBoundaryModels()[patchModelId]->controlMol(*this, td);
}

// ************************************************************************* //

