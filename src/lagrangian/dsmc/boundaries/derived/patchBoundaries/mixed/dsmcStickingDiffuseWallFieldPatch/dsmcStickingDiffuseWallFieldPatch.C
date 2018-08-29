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

#include "dsmcStickingDiffuseWallFieldPatch.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(dsmcStickingDiffuseWallFieldPatch, 0);

addToRunTimeSelectionTable
(
    dsmcPatchBoundary, 
    dsmcStickingDiffuseWallFieldPatch, 
    dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
dsmcStickingDiffuseWallFieldPatch::dsmcStickingDiffuseWallFieldPatch
(
    Time& t,
    const polyMesh& mesh,
    dsmcCloud& cloud,
    const dictionary& dict
)
:
    dsmcPatchBoundary(t, mesh, cloud, dict),
    dsmcFieldPatchBoundary(t, mesh, cloud, dict),
    dsmcStickingDiffuseWallPatch(t, mesh, cloud, dict)
{
    writeInTimeDir_ = false;
    writeInCase_ = false;
    measurePropertiesAtWall_ = true;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

dsmcStickingDiffuseWallFieldPatch::~dsmcStickingDiffuseWallFieldPatch()
{}



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void dsmcStickingDiffuseWallFieldPatch::initialConfiguration()
{}


void dsmcStickingDiffuseWallFieldPatch::calculateProperties()
{}


void dsmcStickingDiffuseWallFieldPatch::controlParticle
(
    dsmcParcel& p,
    dsmcParcel::trackingData& td
)
{
    if(p.isFree())
    {
        const label iD = findIndex(typeIds_, p.typeId());
        
        //- Calculation of the local patch temperature
        const scalar localPatchTemperature = 
            dsmcFieldPatchBoundary::patchLocalTemperature(p);
            
        //- Calculation of the local patch velocity
        const vector& localPatchVelocity = 
            dsmcFieldPatchBoundary::patchLocalVelocity(p);
    
        measurePropertiesBeforeControl(p);
    
        if(iD != -1)
        {
            //- particle considered for adsorption
            const scalar adsorbtionProbability = adsorptionProb(iD);
            
            const label wppIndex = patchId();
            
            const label wppLocalFace = 
                mesh_.boundaryMesh()[wppIndex].whichFace(p.face());
            
            //Info << "Ndot_" << mesh_.time().value() << endl;
    
            if
            (
                adsorbtionProbability > cloud_.rndGen().sample01<scalar>()
             && dsmcStickingWallPatch::isNotSaturated(wppLocalFace)
            )
            {
                //- particle is adsorbed
                dsmcStickingWallPatch::adsorbParticle
                (
                    p,
                    localPatchTemperature
                );
            }
            else 
            {
                //- diffuse reflection
                dsmcDiffuseWallPatch::performDiffuseReflection
                (
                    p,
                    localPatchTemperature,
                    localPatchVelocity
                );
                
                measurePropertiesAfterControl(p);
            }   
        }
        else 
        {
            //- otherwise, it is treated as a diffuse reflection
            dsmcDiffuseWallPatch::performDiffuseReflection
            (
                p,
                localPatchTemperature,
                localPatchVelocity
            );
            
            measurePropertiesAfterControl(p);
        }
    }
    
    //- Separate loop as the particle may have been stuck in the previous loop
    //  If the particle is stuck, consider parcel for release, i.e., desorption
    if(p.isStuck())
    {
        dsmcStickingWallPatch::testForDesorption(p);
    }
    
    //- Update the boundaryMeasurement relative to this sticking patch
    cloud_.boundaryFluxMeasurements().updatenStuckParcelOnPatch
    (
        patchId(),
        dsmcStickingWallPatch::nStuckParcels()
    );
}


void dsmcStickingDiffuseWallFieldPatch::output
(
    const fileName& fixedPathName,
    const fileName& timePath
)
{}


void dsmcStickingDiffuseWallFieldPatch::updateProperties
(
    const dictionary& newDict
)
{
    //- the main properties should be updated first
    updateBoundaryProperties(newDict);
    
    dsmcStickingDiffuseWallPatch::setProperties();
}


} // End namespace Foam

// ************************************************************************* //
