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
    cellMeasurements

Description

\*----------------------------------------------------------------------------*/

#include "cellMeasurements.H"
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
cellMeasurements::cellMeasurements
(
    const polyMesh& mesh,
    dsmcCloud& cloud
)
:
    mesh_(refCast<const fvMesh>(mesh)),
    cloud_(cloud)
{}


//- Construct from mesh, cloud and boolean (dsmcFoam)
cellMeasurements::cellMeasurements
(
    const polyMesh& mesh,
    dsmcCloud& cloud,
    const bool& dummy
)
:
    mesh_(refCast<const fvMesh>(mesh)),
    cloud_(cloud),
    collisionSeparation_(),
    nColls_()
{
    collisionSeparation_.setSize(mesh.nCells(), 0.0);
    nColls_.setSize(mesh.nCells(), 0.0);   
}
// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

cellMeasurements::~cellMeasurements()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void cellMeasurements::clean()
{
    //- clean geometric fields
    
    collisionSeparation_ = 0.0;
    nColls_ = 0.0;
}

void cellMeasurements::reset()
{
    //- reset fields when mesh is edited
    collisionSeparation_.clear();
    collisionSeparation_.setSize(mesh_.nCells(), 0.0);
    
    nColls_.clear();
    nColls_.setSize(mesh_.nCells(), 0.0);
}


void cellMeasurements::updateFields
(
    dsmcParcel& p
)
{
}

// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

}  // End namespace Foam

// ************************************************************************* //
