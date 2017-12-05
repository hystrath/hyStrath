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

#include "dsmcSphericalPatchBoundary.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(dsmcSphericalPatchBoundary, 0);


// * * * * * * * * * * * Protected member functions  * * * * * * * * * * * * //

vector dsmcSphericalPatchBoundary::wallVelocity(const dsmcParcel& p)
{
    const scalar pi = Foam::constant::mathematical::pi;
    
    scalar theta = 0.0; //longitude
    scalar phi = 0.0; //latitude
    scalar radius = 0.0;
    
    radius = mag(centrePoint_ - p.position());
    phi = acos(p.position().y()/radius);
    theta = atan(p.position().x()/p.position().z());
    
    if(p.position().z() < 0 && p.position().x() > 0)
    {
        theta *= -1.0;
        scalar diff = pi/2.0 - theta;
        theta = diff + pi/2.0;
    }
    if(p.position().z() < 0 && p.position().x() < 0)
    {
        theta += pi;
    }
    
//         if(p.position().y() < 500)
//         {
//             Info << "radius = " << radius << endl;
//             Info << "p.position() = " << p.position() << endl;
//             Info << "phi = " << phi*57.2957795 << endl;
//             Info << "theta = " << theta*57.2957795 << endl;
//         }
    
    // add wall velocity
    vector uNew = vector::zero;
    
    scalar linearVelocityXPlane = angularVelocityXPlane_*radius;
    scalar linearVelocityYPlane = angularVelocityYPlane_*radius;
    scalar linearVelocityZPlane = angularVelocityZPlane_*radius;
    
//         if(p.position().y() < 500)
//         {
        uNew.x() = linearVelocityYPlane*cos(theta)*fabs(sin(phi));
        uNew.z() = -linearVelocityYPlane*sin(theta)*fabs(sin(phi));
//         }
//         if( pi < theta <= twoPi)
//         {
//             uNew.x() = -linearVelocityYPlane*cos(theta);
//             uNew.z() = linearVelocityYPlane*sin(theta);
//         }
//         Info << "Ubefore = " << U << endl;
//         if(p.position().y() < 500)
//         {    
//             Info << "uNew = " << uNew << endl;
//         }
    
    return uNew;
}



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
dsmcSphericalPatchBoundary::dsmcSphericalPatchBoundary
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
    angularVelocityXPlane_(propsDict_.lookupOrDefault<scalar>("angularVelocityXPlane", 0.0)),
    angularVelocityYPlane_(propsDict_.lookupOrDefault<scalar>("angularVelocityYPlane", 0.0)),
    angularVelocityZPlane_(propsDict_.lookupOrDefault<scalar>("angularVelocityZPlane", 0.0))
{
    writeInTimeDir_ = false;
    writeInCase_ = false;
    measurePropertiesAtWall_ = true;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

dsmcSphericalPatchBoundary::~dsmcSphericalPatchBoundary()
{}



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void dsmcSphericalPatchBoundary::initialConfiguration()
{}


void dsmcSphericalPatchBoundary::calculateProperties()
{}


void dsmcSphericalPatchBoundary::output
(
    const fileName& fixedPathName,
    const fileName& timePath
)
{}


void dsmcSphericalPatchBoundary::updateProperties(const dictionary& newDict)
{
    //- the main properties should be updated first
    updateBoundaryProperties(newDict);
}


} // End namespace Foam

// ************************************************************************* //
