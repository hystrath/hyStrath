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
    dsmcTimeStepModel

Description

\*----------------------------------------------------------------------------*/

#include "dsmcTimeStepModel.H"
#include "dsmcCloud.H"

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //  
  
    defineTypeNameAndDebug(dsmcTimeStepModel, 0);
    defineRunTimeSelectionTable(dsmcTimeStepModel, fvMesh);
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::dsmcTimeStepModel::initialisenParticles(const scalar value)
{
    forAll(nParticles_, celli)
    {
        nParticles_[celli] = value;
    }
    
    forAll(nParticles_.boundaryField(), patchi)
    {
        forAll(nParticles_.boundaryField()[patchi], facei)
        {
            nParticles_.boundaryFieldRef()[patchi][facei] = value;
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

//- Constructor
Foam::dsmcTimeStepModel::dsmcTimeStepModel
(
    Time& t,
    const polyMesh& mesh,
    dsmcCloud& cloud
)
:
    mesh_(refCast<const fvMesh>(mesh)),
    nParticles_
    (
        IOobject
        (
            "nParticles",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("nParticles", dimless, 1.0)
    )
{}


// * * * * * * * * * * * * * * * * Selector  * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::dsmcTimeStepModel>
Foam::dsmcTimeStepModel::New
(
    Time& t,
    const polyMesh& mesh,
    dsmcCloud& cloud
)
{
    word timeStepModel =
        cloud.particleProperties().lookupOrDefault<word>
        (
            "timeStepModel",
            "constant"
        );
      
    timeStepModel = "dsmc" + static_cast<word>(std::toupper(timeStepModel[0]))
        + timeStepModel.substr(1) + "TimeStepModel";

    Info<< "Selecting the time-step model:" << tab << timeStepModel 
        << "\n" << endl;

    fvMeshConstructorTable::iterator cstrIter =
        fvMeshConstructorTablePtr_->find(timeStepModel);

    if (cstrIter == fvMeshConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "dsmcTimeStepModel::New(Time&, const polyMesh&, dsmcCloud&)"
        )   << "Unknown time-step model type "
            << timeStepModel << endl << endl
            << "Valid time-step model types are : " << endl
            << fvMeshConstructorTablePtr_->toc()
            << exit(FatalError);
    }

    return autoPtr<dsmcTimeStepModel>(cstrIter()(t, mesh, cloud));
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::dsmcTimeStepModel::~dsmcTimeStepModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //





// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
