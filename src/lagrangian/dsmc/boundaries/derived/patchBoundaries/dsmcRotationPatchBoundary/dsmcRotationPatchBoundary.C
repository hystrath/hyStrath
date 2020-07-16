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

#include "dsmcRotationPatchBoundary.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(dsmcRotationPatchBoundary, 0);


// * * * * * * * * * * * Protected member functions  * * * * * * * * * * * * //

vector dsmcRotationPatchBoundary::wallVelocity(const dsmcParcel& p)
{
    //- Calculation of the wall velocity uNew to be added to U
    vector uNew =
        (
          (p.position() - centrePoint_)/mag((p.position() - centrePoint_))
        ) ^ rotationAxis_;

    uNew /= mag(uNew);

    uNew *= wallVelocity_;

    return uNew;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
dsmcRotationPatchBoundary::dsmcRotationPatchBoundary
(
    Time& t,
    const polyMesh& mesh,
    dsmcCloud& cloud,
    const dictionary& dict
)
:
    dsmcPatchBoundary(t, mesh, cloud, dict),
    propsDict_(dict.subDict(typeName + "Properties")),
    centrePoint_(propsDict_.lookupOrDefault<vector>("centrePoint", vector::zero)),
    rotationAxis_(propsDict_.lookupOrDefault<vector>("rotationAxis", vector(1, 0, 0))),
    wallVelocity_(propsDict_.lookupOrDefault<scalar>("wallVelocity", 0.0))
{
    writeInTimeDir_ = false;
    writeInCase_ = false;
    measurePropertiesAtWall_ = true;

    rotationAxis_ /= mag(rotationAxis_);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

dsmcRotationPatchBoundary::~dsmcRotationPatchBoundary()
{}



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void dsmcRotationPatchBoundary::initialConfiguration()
{}


void dsmcRotationPatchBoundary::calculateProperties()
{}


void dsmcRotationPatchBoundary::output
(
    const fileName& fixedPathName,
    const fileName& timePath
)
{}


void dsmcRotationPatchBoundary::updateProperties(const dictionary& newDict)
{
    //- the main properties should be updated first
    updateBoundaryProperties(newDict);
}


} // End namespace Foam

// ************************************************************************* //
