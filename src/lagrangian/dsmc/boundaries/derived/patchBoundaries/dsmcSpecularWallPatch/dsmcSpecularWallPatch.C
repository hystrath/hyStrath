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

#include "dsmcSpecularWallPatch.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(dsmcSpecularWallPatch, 0);

addToRunTimeSelectionTable
(
    dsmcPatchBoundary,
    dsmcSpecularWallPatch,
    dictionary
);

// * * * * * * * * * * * Protected member functions  * * * * * * * * * * * * //

void dsmcSpecularWallPatch::performSpecularReflection(dsmcParcel& p)
{
    vector& U = p.U();

    vector nw = p.normal();
    nw /= mag(nw);

    const scalar& U_dot_nw = U & nw;

    if (U_dot_nw > 0.0)
    {
        U -= 2.0*U_dot_nw*nw;
    }

    cloud_.porousMeas().specularInteraction(p, nw);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
dsmcSpecularWallPatch::dsmcSpecularWallPatch
(
    Time& t,
    const polyMesh& mesh,
    dsmcCloud& cloud,
    const dictionary& dict
)
:
    dsmcPatchBoundary(t, mesh, cloud, dict),
    propsDict_(dict.subDict(typeName + "Properties"))
{
    writeInTimeDir_ = false;
    writeInCase_ = false;
    measurePropertiesAtWall_ = true;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

dsmcSpecularWallPatch::~dsmcSpecularWallPatch()
{}



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
void dsmcSpecularWallPatch::initialConfiguration()
{}


void dsmcSpecularWallPatch::calculateProperties()
{}


void dsmcSpecularWallPatch::controlParticle
(
    dsmcParcel& p,
    dsmcParcel::trackingData& td
)
{
    measurePropertiesBeforeControl(p);

    performSpecularReflection(p);

    measurePropertiesAfterControl(p);
}


void dsmcSpecularWallPatch::output
(
    const fileName& fixedPathName,
    const fileName& timePath
)
{}


void dsmcSpecularWallPatch::updateProperties(const dictionary& newDict)
{
    //- the main properties should be updated first
    updateBoundaryProperties(newDict);
}


} // End namespace Foam

// ************************************************************************* //
