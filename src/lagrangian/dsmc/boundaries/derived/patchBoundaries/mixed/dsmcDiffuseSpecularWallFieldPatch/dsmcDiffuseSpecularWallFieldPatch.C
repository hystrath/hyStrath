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

#include "dsmcDiffuseSpecularWallFieldPatch.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(dsmcDiffuseSpecularWallFieldPatch, 0);

addToRunTimeSelectionTable
(
    dsmcPatchBoundary,
    dsmcDiffuseSpecularWallFieldPatch,
    dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
dsmcDiffuseSpecularWallFieldPatch::dsmcDiffuseSpecularWallFieldPatch
(
    Time& t,
    const polyMesh& mesh,
    dsmcCloud& cloud,
    const dictionary& dict
)
:
    dsmcPatchBoundary(t, mesh, cloud, dict),
    dsmcFieldPatchBoundary(t, mesh, cloud, dict),
    dsmcDiffuseSpecularWallPatch(t, mesh, cloud, dict),
    propsDict_(dict.subDict(typeName + "Properties"))
{
    writeInTimeDir_ = false;
    writeInCase_ = false;
    measurePropertiesAtWall_ = true;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

dsmcDiffuseSpecularWallFieldPatch::~dsmcDiffuseSpecularWallFieldPatch()
{}



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void dsmcDiffuseSpecularWallFieldPatch::initialConfiguration()
{}


void dsmcDiffuseSpecularWallFieldPatch::calculateProperties()
{}


void dsmcDiffuseSpecularWallFieldPatch::controlParticle
(
    dsmcParcel& p,
    dsmcParcel::trackingData& td
)
{
    measurePropertiesBeforeControl(p);

    if (diffuseFraction_ > cloud_.rndGen().sample01<scalar>())
    {
        //- Calculation of the local patch temperature
        const scalar& localPatchTemperature =
            dsmcFieldPatchBoundary::patchLocalTemperature(p);

        //- Calculation of the local patch velocity
        const vector& localPatchVelocity =
            dsmcFieldPatchBoundary::patchLocalVelocity(p);

        //- Diffuse reflection with local info on T and U
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


void dsmcDiffuseSpecularWallFieldPatch::output
(
    const fileName& fixedPathName,
    const fileName& timePath
)
{}


void dsmcDiffuseSpecularWallFieldPatch::updateProperties
(
    const dictionary& newDict
)
{
    //- the main properties should be updated first
    updateBoundaryProperties(newDict);
}


} // End namespace Foam

// ************************************************************************* //
