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
    dsmcVariableTimeStepModel

Description

\*----------------------------------------------------------------------------*/

#include "addToRunTimeSelectionTable.H"
#include "dsmcVariableTimeStepModel.H"
#include "dsmcCloud.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(dsmcVariableTimeStepModel, 0);

    addToRunTimeSelectionTable
    (
        dsmcTimeStepModel, 
        dsmcVariableTimeStepModel,
        fvMesh
    );

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void dsmcVariableTimeStepModel::findRefCell()
{
    // Find the cell with minimum volume
    // This cell will serve as a reference cell to set the 
    // nParticle/timeStep ratio
    
    const scalarField& volumeCells = mesh_.V();
    scalar minVolume = gMin(volumeCells);
    
    if (Pstream::parRun())
    {
        reduce(minVolume, minOp<scalar>());
    }

    forAll(nParticles_, celli)
    {
        if (mag(volumeCells[celli] - minVolume) < SMALL)
        {
            refCell_ = celli;
            break;
        }
    }
}


void dsmcVariableTimeStepModel::updatenParticles()
{
    findRefCell();
    
    const scalarField& volumeCells = mesh_.V();
    const scalar minVolume = volumeCells[refCell_];
    
    const scalar nParticleRef = nParticles_[refCell_];
        
    forAll(nParticles_, celli)
    {
        nParticles_[celli] = nParticleRef*volumeCells[celli]/minVolume;
    }
    
    forAll(nParticles_.boundaryField(), patchi)
    {
        fvPatchScalarField& pnParticles = 
            nParticles_.boundaryFieldRef()[patchi];
        
        forAll(pnParticles, facei)
        {
            pnParticles[facei] = 
                nParticles_[mesh_.boundaryMesh()[patchi].faceCells()[facei]];
        }
    }
}


void dsmcVariableTimeStepModel::updateTimeStep()
{
    const scalar nParticleTimeStepRatio = 
        nParticles_[refCell_]/deltaT_[refCell_];
    
    forAll(deltaT_, celli)
    {
        deltaT_[celli] = nParticles_[celli]/nParticleTimeStepRatio;
    }
}


void dsmcVariableTimeStepModel::updateVariableTimeStepMethod()
{
    updatenParticles();
    
    updateTimeStep();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

//- Constructor
dsmcVariableTimeStepModel::dsmcVariableTimeStepModel
(
    Time& t,
    const polyMesh& mesh,
    dsmcCloud& cloud
)
:
    dsmcTimeStepModel(t, mesh, cloud),
    cloud_(cloud),
    refCell_(-1),
    deltaT_
    (
        IOobject
        (
            "deltaT",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("deltaT", dimTime, deltaTValueOrg())
    )
{
    nParticles_.writeOpt() = IOobject::AUTO_WRITE;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

dsmcVariableTimeStepModel::~dsmcVariableTimeStepModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void dsmcVariableTimeStepModel::checkTimeStepModelInputs()
{
    const bool nParticlesFromFile = 
        cloud_.particleProperties().lookupOrDefault<bool>
        (
            "nEquivalentParticlesFromFile",
            false
        );
    
    if(nParticlesFromFile)
    {
        nParticles_.regIOobject::read();
        
        findRefCell();
    }
    else
    {
        updatenParticles();
    }
    
    updateTimeStep();
    
    writeTimeStepModelInfo();
}


void dsmcVariableTimeStepModel::update()
{
    updateVariableTimeStepMethod();
}


void dsmcVariableTimeStepModel::writeTimeStepModelInfo() const
{
    Info<< "Variable time-step model:" << nl
        << "- Initial time-step [sec]" << tab
        << dsmcTimeStepModel::deltaTValue(0) << nl
        << endl;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
