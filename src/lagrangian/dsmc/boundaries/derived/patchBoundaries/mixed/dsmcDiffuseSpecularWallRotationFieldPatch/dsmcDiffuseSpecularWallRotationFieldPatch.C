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

Description

\*---------------------------------------------------------------------------*/

#include "dsmcDiffuseSpecularWallRotationFieldPatch.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(dsmcDiffuseSpecularWallRotationFieldPatch, 0);

addToRunTimeSelectionTable
(
    dsmcPatchBoundary,
    dsmcDiffuseSpecularWallRotationFieldPatch,
    dictionary
);


// * * * * * * * * * * * Protected member functions  * * * * * * * * * * * * //

void dsmcDiffuseSpecularWallRotationFieldPatch::setProperties()
{
    dsmcDiffuseSpecularWallPatch::setProperties();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
dsmcDiffuseSpecularWallRotationFieldPatch::dsmcDiffuseSpecularWallRotationFieldPatch
(
    Time& t,
    const polyMesh& mesh,
    dsmcCloud& cloud,
    const dictionary& dict
)
:
    dsmcPatchBoundary(t, mesh, cloud, dict),
    dsmcRotationPatchBoundary(t, mesh, cloud, dict),
    dsmcFieldPatchBoundary(t, mesh, cloud, dict),
    dsmcDiffuseSpecularWallPatch(t, mesh, cloud, dict),
    propsDict_(dict.subDict(typeName + "Properties"))

{
    writeInTimeDir_ = false;
    writeInCase_ = false;
    measurePropertiesAtWall_ = true;

    setProperties();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

dsmcDiffuseSpecularWallRotationFieldPatch::~dsmcDiffuseSpecularWallRotationFieldPatch()
{}



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void dsmcDiffuseSpecularWallRotationFieldPatch::initialConfiguration()
{}


void dsmcDiffuseSpecularWallRotationFieldPatch::calculateProperties()
{}


void dsmcDiffuseSpecularWallRotationFieldPatch::controlParticle
(
    dsmcParcel& p,
    dsmcParcel::trackingData& td
)
{
    measurePropertiesBeforeControl(p);

    if (diffuseFraction() > cloud_.rndGen().sample01<scalar>())
    {
        //- Calculation of the wall velocity to be added to U
        const vector& localPatchVelocity =
            dsmcRotationPatchBoundary::wallVelocity(p);

        //- Calculation of the local patch temperature
        const scalar& localPatchTemperature =
            dsmcFieldPatchBoundary::patchLocalTemperature(p);

        //- Diffuse reflection
        dsmcDiffuseWallPatch::performDiffuseReflection
        (
            p,
            localPatchTemperature,
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


void dsmcDiffuseSpecularWallRotationFieldPatch::output
(
    const fileName& fixedPathName,
    const fileName& timePath
)
{}


void dsmcDiffuseSpecularWallRotationFieldPatch::updateProperties(const dictionary& newDict)
{
    //- the main properties should be updated first
    updateBoundaryProperties(newDict);

    setProperties();
}


} // End namespace Foam

// ************************************************************************* //
