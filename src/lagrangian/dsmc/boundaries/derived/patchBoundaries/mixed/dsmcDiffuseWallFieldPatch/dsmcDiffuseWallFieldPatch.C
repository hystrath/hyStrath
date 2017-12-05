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

#include "dsmcDiffuseWallFieldPatch.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(dsmcDiffuseWallFieldPatch, 0);

addToRunTimeSelectionTable
(
    dsmcPatchBoundary, 
    dsmcDiffuseWallFieldPatch, 
    dictionary
);


// * * * * * * * * * * * Protected member functions  * * * * * * * * * * * * //

void dsmcDiffuseWallFieldPatch::performDiffuseFieldReflection(dsmcParcel& p)
{
    //- Calculation of the wall velocity to be added to U
    const vector& localPatchVelocity = 
        dsmcFieldPatchBoundary::patchLocalVelocity(p);
    
    //- Calculation of the local patch temperature
    const scalar& localPatchTemperature = 
        dsmcFieldPatchBoundary::patchLocalTemperature(p);
    
    //- Diffuse reflection with local info on T and U
    dsmcDiffuseWallPatch::performDiffuseReflection
    (
        p, 
        localPatchTemperature, 
        localPatchVelocity
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
dsmcDiffuseWallFieldPatch::dsmcDiffuseWallFieldPatch
(
    Time& t,
    const polyMesh& mesh,
    dsmcCloud& cloud,
    const dictionary& dict
)
:
    dsmcPatchBoundary(t, mesh, cloud, dict),
    dsmcFieldPatchBoundary(t, mesh, cloud, dict),
    dsmcDiffuseWallPatch(t, mesh, cloud, dict),
    propsDict_(dict.subDict(dsmcDiffuseWallPatch::typeName + "Properties"))
{
    writeInTimeDir_ = false;
    writeInCase_ = false;
    measurePropertiesAtWall_ = true;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

dsmcDiffuseWallFieldPatch::~dsmcDiffuseWallFieldPatch()
{}



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void dsmcDiffuseWallFieldPatch::initialConfiguration()
{}


void dsmcDiffuseWallFieldPatch::calculateProperties()
{}


void dsmcDiffuseWallFieldPatch::controlParticle
(
    dsmcParcel& p, 
    dsmcParcel::trackingData& td
)
{
    measurePropertiesBeforeControl(p);

    performDiffuseFieldReflection(p);

    measurePropertiesAfterControl(p);
}


void dsmcDiffuseWallFieldPatch::output
(
    const fileName& fixedPathName,
    const fileName& timePath
)
{}


void dsmcDiffuseWallFieldPatch::updateProperties(const dictionary& newDict)
{
    //- the main properties should be updated first
    updateBoundaryProperties(newDict);
}


} // End namespace Foam

// ************************************************************************* //
