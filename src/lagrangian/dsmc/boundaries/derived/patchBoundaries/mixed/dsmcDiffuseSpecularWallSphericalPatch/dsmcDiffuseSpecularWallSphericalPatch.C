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

#include "dsmcDiffuseSpecularWallSphericalPatch.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(dsmcDiffuseSpecularWallSphericalPatch, 0);

addToRunTimeSelectionTable
(
    dsmcPatchBoundary, 
    dsmcDiffuseSpecularWallSphericalPatch, 
    dictionary
);


// * * * * * * * * * * * Protected member functions  * * * * * * * * * * * * //

void dsmcDiffuseSpecularWallSphericalPatch::setProperties()
{
    dsmcDiffuseSpecularWallPatch::setProperties();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
dsmcDiffuseSpecularWallSphericalPatch::dsmcDiffuseSpecularWallSphericalPatch
(
    Time& t,
    const polyMesh& mesh,
    dsmcCloud& cloud,
    const dictionary& dict
)
:
    dsmcPatchBoundary(t, mesh, cloud, dict),
    dsmcSphericalPatchBoundary(t, mesh, cloud, dict),
    dsmcDiffuseSpecularWallPatch(t, mesh, cloud, dict),
    propsDict_(dict.subDict(typeName + "Properties"))

{
    writeInTimeDir_ = false;
    writeInCase_ = false;
    measurePropertiesAtWall_ = true;
    
    setProperties();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

dsmcDiffuseSpecularWallSphericalPatch::~dsmcDiffuseSpecularWallSphericalPatch()
{}



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void dsmcDiffuseSpecularWallSphericalPatch::initialConfiguration()
{}


void dsmcDiffuseSpecularWallSphericalPatch::calculateProperties()
{}


void dsmcDiffuseSpecularWallSphericalPatch::controlParticle(dsmcParcel& p, dsmcParcel::trackingData& td)
{
    measurePropertiesBeforeControl(p);

    //- Calculation of the wall velocity to be added to U
    const vector& localPatchVelocity = 
        dsmcSphericalPatchBoundary::wallVelocity(p);

    if (diffuseFraction() > cloud_.rndGen().scalar01())
    {
        //- Diffuse reflection
        dsmcDiffuseWallPatch::performDiffuseReflection
        (
            p,
            0,
            localPatchVelocity
        );
    }
    else
    {
        //- Specular reflection
        dsmcSpecularWallPatch::performSpecularReflection(p);
    }

    measurePropertiesAfterControl(p);
}

void dsmcDiffuseSpecularWallSphericalPatch::output
(
    const fileName& fixedPathName,
    const fileName& timePath
)
{}


void dsmcDiffuseSpecularWallSphericalPatch::updateProperties(const dictionary& newDict)
{
    //- the main properties should be updated first
    updateBoundaryProperties(newDict);
}


} // End namespace Foam

// ************************************************************************* //
