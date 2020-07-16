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


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
