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

#include "dsmcDiffuseSpecularWallPatch.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(dsmcDiffuseSpecularWallPatch, 0);

addToRunTimeSelectionTable
(
    dsmcPatchBoundary, 
    dsmcDiffuseSpecularWallPatch, 
    dictionary
);


// * * * * * * * * * * * Protected member functions  * * * * * * * * * * * * //

void dsmcDiffuseSpecularWallPatch::setProperties()
{
    dsmcDiffuseWallPatch::setProperties();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
dsmcDiffuseSpecularWallPatch::dsmcDiffuseSpecularWallPatch
(
    Time& t,
    const polyMesh& mesh,
    dsmcCloud& cloud,
    const dictionary& dict
)
:
    dsmcPatchBoundary(t, mesh, cloud, dict),
    dsmcSpecularWallPatch(t, mesh, cloud, dict),
    dsmcDiffuseWallPatch(t, mesh, cloud, dict),
    propsDict_(dict.subDict(typeName + "Properties")),
    diffuseFraction_(readScalar(propsDict_.lookup("diffuseFraction")))
{
    writeInTimeDir_ = false;
    writeInCase_ = false;
    measurePropertiesAtWall_ = true;

    setProperties();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

dsmcDiffuseSpecularWallPatch::~dsmcDiffuseSpecularWallPatch()
{}



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
void dsmcDiffuseSpecularWallPatch::initialConfiguration()
{}


void dsmcDiffuseSpecularWallPatch::calculateProperties()
{}


void dsmcDiffuseSpecularWallPatch::controlParticle
(
    dsmcParcel& p, 
    dsmcParcel::trackingData& td
)
{
    measurePropertiesBeforeControl(p);

    if (diffuseFraction_ > cloud_.rndGen().sample01<scalar>())
    {
        //- Diffuse reflection
        dsmcDiffuseWallPatch::performDiffuseReflection(p);
    }
    else
    {
        //- Specular reflection
        dsmcSpecularWallPatch::performSpecularReflection(p);
    }

    measurePropertiesAfterControl(p);
}


void dsmcDiffuseSpecularWallPatch::output
(
    const fileName& fixedPathName,
    const fileName& timePath
)
{}


void dsmcDiffuseSpecularWallPatch::updateProperties(const dictionary& newDict)
{
    //- the main properties should be updated first
    updateBoundaryProperties(newDict);

    setProperties();
}


} // End namespace Foam

// ************************************************************************* //
