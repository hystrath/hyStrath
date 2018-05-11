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
    dsmcInstantaneousPorousMediumMeasurements

Description

\*----------------------------------------------------------------------------*/

#include "addToRunTimeSelectionTable.H"
#include "dsmcInstantaneousPorousMediumMeasurements.H"
#include "dsmcCloud.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(dsmcInstantaneousPorousMediumMeasurements, 0);

    addToRunTimeSelectionTable
    (
        porousMeasurements, 
        dsmcInstantaneousPorousMediumMeasurements,
        fvMesh
    );

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void dsmcInstantaneousPorousMediumMeasurements::writeInstantaneousPorousMeasurementsInfo() const
{
    Info<< nl << "Porous media measurements - instantaneous:" << nl
        << "- tracking probability" << tab << trackingProbability_ << nl
        << "- dimensionality" << tab << dimensionality_ << nl
        << endl;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

//- Constructor
dsmcInstantaneousPorousMediumMeasurements::dsmcInstantaneousPorousMediumMeasurements
(
    Time& t,
    const polyMesh& mesh,
    dsmcCloud& cloud
)
:
    porousMeasurements(t, mesh, cloud),
    mesh_(refCast<const fvMesh>(mesh)),
    cloud_(cloud),
    trackingProbability_(0),
    dimensionality_(3),
    initialSimulationTime_(mesh_.time().value()),
    mediumTotalDistanceTravelled_(0.0)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

dsmcInstantaneousPorousMediumMeasurements::~dsmcInstantaneousPorousMediumMeasurements()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void dsmcInstantaneousPorousMediumMeasurements::checkPorousMeasurementsInputs()
{
    if(cloud_.particleProperties().isDict("tracerProperties"))
    {
        trackingProbability_ = cloud_.particleProperties()
            .subDict("tracerProperties")
            .lookupOrDefault<scalar>("trackingProbability", 0.1);
            
        dimensionality_ = cloud_.particleProperties()
            .subDict("tracerProperties")
            .lookupOrDefault<label>("dimensionality", 3);
    }
    else
    {
        WarningIn
        (
            "dsmcInstantaneousPorousMediumMeasurements::"
            "checkPorousMeasurementsInputs()"
        )
        << "Subdictionary tracerProperties is missing in "
           "constant/dsmcProperties"
        << endl;
        
        trackingProbability_ = 0.1;
        dimensionality_ = 3;
    }
    
    forAllIter(dsmcCloud, cloud_, iter)
    {
        if(not iter().isTracked())
        {
            if(trackingProbability_ > cloud_.rndGen().scalar01())
            {
                iter().setTracked
                (
                    true,
                    -1,
                    mesh_.time().value(), 
                    iter().position()
                );
            }
        }
    }
        
    writeInstantaneousPorousMeasurementsInfo();
}


void dsmcInstantaneousPorousMediumMeasurements::update()
{}


void dsmcInstantaneousPorousMediumMeasurements::specularInteraction
(
    dsmcParcel& p,
    const vector& nw
)
{
    if (p.isTracked())
    {
        if (p.tracked().inPatchId() == -1)
        {
            p.tracked().updateDistanceTravelled(p.position());
        
            p.tracked().performSpecularReflectionOnDistanceTravelled(nw);
        }
    } 
}


void dsmcInstantaneousPorousMediumMeasurements::diffuseInteraction
(
    dsmcParcel& p
)
{}


void dsmcInstantaneousPorousMediumMeasurements::cyclicMembraneInteraction
(
    dsmcParcel& p,
    const vector& orgPos
)
{
    if (p.isTracked())
    {
        if (p.tracked().inPatchId() == -1)
        {
            p.tracked().updateDistanceTravelled(orgPos);
            
            p.tracked().updateCurrentPosition(p.position());
        }
    }
}

    
void
dsmcInstantaneousPorousMediumMeasurements::writePorousMeasurementsInfo() const
{
    vector iMSDvector = vector::zero;
    label cmpt = 0;

    forAllIter(dsmcCloud, cloud_, iter)
    {
        dsmcParcel& p = iter();
        
        if (p.isTracked())
        {
            if (p.tracked().inPatchId() == -1)
            {
                p.tracked().updateDistanceTravelled(p.position());
                
                p.tracked().updateMeanSquaredDisplacement();
                
                iMSDvector += p.tracked().meanSquaredDisplacementVector();
                
                cmpt++;
            }
        }
    }

    reduce(cmpt, sumOp<label>());
    reduce(iMSDvector, sumOp<vector>());
    
    const scalar iMSD = iMSDvector[0] + iMSDvector[1] + iMSDvector[2];
    
    const scalar nTracers = cmpt + SMALL;
    
    const scalar deltaSimulationTime = mesh_.time().value() 
        - initialSimulationTime_;
      
    Info<< "    Mean squared displacement       = " 
        << iMSD/nTracers << " (" << cmpt << ")" << nl
        << "    Mean squared displacement_x     = " 
        << iMSDvector[0]/nTracers << nl
        << "    Mean squared displacement_y     = " 
        << iMSDvector[1]/nTracers << nl
        << "    Mean squared displacement_z     = " 
        << iMSDvector[2]/nTracers << nl
        << "    Effective diffusivity (approx.) = "
        << iMSD/(2.0*dimensionality_*deltaSimulationTime*nTracers)
        << endl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
