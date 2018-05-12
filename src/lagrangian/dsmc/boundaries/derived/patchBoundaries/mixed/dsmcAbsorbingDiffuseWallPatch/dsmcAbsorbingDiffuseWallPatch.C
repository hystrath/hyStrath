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

#include "dsmcAbsorbingDiffuseWallPatch.H"
#include "addToRunTimeSelectionTable.H"
#include "fvc.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(dsmcAbsorbingDiffuseWallPatch, 0);

addToRunTimeSelectionTable
(
    dsmcPatchBoundary, 
    dsmcAbsorbingDiffuseWallPatch, 
    dictionary
);


// * * * * * * * * * * * Protected member functions  * * * * * * * * * * * * //

void dsmcAbsorbingDiffuseWallPatch::setProperties()
{
    dsmcDiffuseWallPatch::setProperties();
    dsmcAbsorbingWallPatch::setProperties();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
dsmcAbsorbingDiffuseWallPatch::dsmcAbsorbingDiffuseWallPatch
(
    Time& t,
    const polyMesh& mesh,
    dsmcCloud& cloud,
    const dictionary& dict
)
:
    dsmcPatchBoundary(t, mesh, cloud, dict),
    dsmcDiffuseWallPatch(t, mesh, cloud, dict),
    dsmcAbsorbingWallPatch(t, mesh, cloud, dict)
{
    writeInTimeDir_ = false;
    writeInCase_ = false;
    measurePropertiesAtWall_ = true;
    
    dsmcDiffuseWallPatch::setProperties();
    dsmcAbsorbingWallPatch::setProperties();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

dsmcAbsorbingDiffuseWallPatch::~dsmcAbsorbingDiffuseWallPatch()
{}



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void dsmcAbsorbingDiffuseWallPatch::initialConfiguration()
{}


void dsmcAbsorbingDiffuseWallPatch::calculateProperties()
{}


void dsmcAbsorbingDiffuseWallPatch::controlParticle
(
    dsmcParcel& p, 
    dsmcParcel::trackingData& td
)
{
    measurePropertiesBeforeControl(p);
    
    const label iD = findIndex(typeIds_, p.typeId());
    
    if(iD != -1) 
    {
        //- particle considered for absorption
        const scalar absorptionProbability = absorptionProb(iD);
        
        const label wppIndex = patchId();
        
        const label wppLocalFace = 
            mesh_.boundaryMesh()[wppIndex].whichFace(p.face());
        
        if
        (
            absorptionProbability > cloud_.rndGen().sample01<scalar>()
         && dsmcAbsorbingWallPatch::isNotSaturated(wppIndex, wppLocalFace)
        )
        {
            //- absorb particle
            absorbParticle(wppIndex, wppLocalFace, td);
        }
        else
        {
            //- diffuse reflection
            dsmcDiffuseWallPatch::performDiffuseReflection(p);
            
            measurePropertiesAfterControl(p);
        }   
    }
    else
    {
        //- otherwise, it is treated as a diffuse reflection
        dsmcDiffuseWallPatch::performDiffuseReflection(p);
        
        measurePropertiesAfterControl(p);
    }
}


void dsmcAbsorbingDiffuseWallPatch::output
(
    const fileName& fixedPathName,
    const fileName& timePath
)
{}


void dsmcAbsorbingDiffuseWallPatch::updateProperties(const dictionary& newDict)
{
    //- the main properties should be updated first
    updateBoundaryProperties(newDict);

    setProperties();
}


} // End namespace Foam

// ************************************************************************* //
