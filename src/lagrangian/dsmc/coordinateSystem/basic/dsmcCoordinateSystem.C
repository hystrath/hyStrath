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
    dsmcCoordinateSystem

Description

\*----------------------------------------------------------------------------*/

#include "dsmcCoordinateSystem.H"
#include "dsmcCloud.H"

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //  
  
    defineTypeNameAndDebug(dsmcCoordinateSystem, 0);
    defineRunTimeSelectionTable(dsmcCoordinateSystem, fvMesh);
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

//- Constructor
Foam::dsmcCoordinateSystem::dsmcCoordinateSystem
(
    Time& t,
    const polyMesh& mesh,
    dsmcCloud& cloud
)
:
    mesh_(refCast<const fvMesh>(mesh)),
    nParticlesOrg_
    (
        readScalar(cloud.particleProperties().lookup("nEquivalentParticles"))
    ),
    timeStepModel_(dsmcTimeStepModel::New(t, mesh, cloud))
{
    dtModel().initialisenParticles(nParticlesOrg_);
    
    dtModel().checkTimeStepModelInputs();
}


// * * * * * * * * * * * * * * * * Selector  * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::dsmcCoordinateSystem>
Foam::dsmcCoordinateSystem::New
(
    Time& t,
    const polyMesh& mesh,
    dsmcCloud& cloud
)
{
   const word& coordSystem =
      cloud.particleProperties().lookupOrDefault<word>
      (
          "coordinateSystem",
          "dsmcCartesian"
      );
    
    Info<< "Selecting the coordinate system model:" << tab << coordSystem 
        << "\n" << endl;

    fvMeshConstructorTable::iterator cstrIter =
        fvMeshConstructorTablePtr_->find(coordSystem);

    if (cstrIter == fvMeshConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "dsmcCoordinateSystem::New(Time&, const polyMesh& mesh,"
            "dsmcCloud& cloud)"
        )   << "Unknown coordinate system type "
            << coordSystem << endl << endl
            << "Valid coordinate system types are : " << endl
            << fvMeshConstructorTablePtr_->toc()
            << exit(FatalError);
    }

    return autoPtr<dsmcCoordinateSystem>(cstrIter()(t, mesh, cloud));
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::dsmcCoordinateSystem::~dsmcCoordinateSystem()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::dsmcCoordinateSystem::recalculatepRWF
(
    const label patchI,
    const label faceI
) const
{
    return 1.0;    
}


Foam::scalar Foam::dsmcCoordinateSystem::recalculateRWF
(
    const label cellI, 
    const bool mixedRWMethod
) const
{
    return 1.0;    
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
