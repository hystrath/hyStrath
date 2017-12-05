/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2007 OpenCFD Ltd.
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

Description

\*---------------------------------------------------------------------------*/

#include "dsmcFieldPatchBoundary.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(dsmcFieldPatchBoundary, 0);


// * * * * * * * * * * * Protected member functions  * * * * * * * * * * * * //

void dsmcFieldPatchBoundary::readPatchFields()
{
    //- Temperature field
    tmp<volScalarField> tboundaryT
    (
        new volScalarField
        (
            IOobject
            (
                "boundaryT",
                mesh_.time().timeName(),
                mesh_,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            ),
            mesh_
        )
    );
    
    volScalarField& boundaryT = tboundaryT.ref();
    
    boundaryT_ = boundaryT.boundaryField()[patchId()];

    cloud_.boundaryFluxMeasurements().setBoundaryT(patchId(), boundaryT_);
    
    //- Velocity field
    tmp<volVectorField> tboundaryU
    (
        new volVectorField
        (
            IOobject
            (
                "boundaryU",
                mesh_.time().timeName(),
                mesh_,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            ),
            mesh_
        )
    );
    
    volVectorField& boundaryU = tboundaryU.ref();
    
    boundaryU_ = boundaryU.boundaryField()[patchId()];
    
    cloud_.boundaryFluxMeasurements().setBoundaryU(patchId(), boundaryU_);
}

/*scalar dsmcFieldPatchBoundary::patchLocalTemperature(const dsmcParcel& p)
{
    const label wppIndex = p.patch(p.face());

    const polyPatch& patch = mesh_.boundaryMesh()[wppIndex];

    const label wppLocalFace = patch.whichFace(p.face());
    
    return boundaryT_.boundaryField()[wppIndex][wppLocalFace];
}


vector dsmcFieldPatchBoundary::patchLocalVelocity(const dsmcParcel& p)
{
    const label wppIndex = p.patch(p.face());

    const polyPatch& patch = mesh_.boundaryMesh()[wppIndex];

    const label wppLocalFace = patch.whichFace(p.face());
    
    return boundaryU_.boundaryField()[wppIndex][wppLocalFace];
}*/


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
dsmcFieldPatchBoundary::dsmcFieldPatchBoundary
(
    Time& t,
    const polyMesh& mesh,
    dsmcCloud& cloud,
    const dictionary& dict
)
:
    dsmcPatchBoundary(t, mesh, cloud, dict),
    boundaryT_(mesh_.boundaryMesh()[patchId()].size(), 0.0),
    boundaryU_(mesh_.boundaryMesh()[patchId()].size(), vector::zero)
    /*boundaryT_
    (
        volScalarField
        (
            IOobject
            (
                "boundaryT",
                mesh_.time().timeName(),
                mesh_,
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            mesh_
        )
    ),
    boundaryU_
    (
        volVectorField
        (
            IOobject
            (
                "boundaryU",
                mesh_.time().timeName(),
                mesh_,
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            mesh_
        )
    )*/
{
    writeInTimeDir_ = false;
    writeInCase_ = false;
    measurePropertiesAtWall_ = true;
    
    readPatchFields();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

dsmcFieldPatchBoundary::~dsmcFieldPatchBoundary()
{}



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void dsmcFieldPatchBoundary::initialConfiguration()
{}


void dsmcFieldPatchBoundary::calculateProperties()
{}


void dsmcFieldPatchBoundary::output
(
    const fileName& fixedPathName,
    const fileName& timePath
)
{}


void dsmcFieldPatchBoundary::updateProperties(const dictionary& newDict)
{
    //- the main properties should be updated first
    updateBoundaryProperties(newDict);
    
    readPatchFields();
}


//- Access

scalar dsmcFieldPatchBoundary::patchLocalTemperature(const dsmcParcel& p) const
{
    const label wppLocalFace = 
        mesh_.boundaryMesh()[patchId()].whichFace(p.face());
    
    return boundaryT_[wppLocalFace];
}


const vector&
dsmcFieldPatchBoundary::patchLocalVelocity(const dsmcParcel& p) const
{
    const label wppLocalFace = 
        mesh_.boundaryMesh()[patchId()].whichFace(p.face());
    
    return boundaryU_[wppLocalFace];
}


} // End namespace Foam

// ************************************************************************* //
