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

#include "dsmcMixedDiffuseSpecularWallRotationPatch.H"
#include "addToRunTimeSelectionTable.H"
#include "fvc.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(dsmcMixedDiffuseSpecularWallRotationPatch, 0);

addToRunTimeSelectionTable(dsmcPatchBoundary, dsmcMixedDiffuseSpecularWallRotationPatch, dictionary);



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
dsmcMixedDiffuseSpecularWallRotationPatch::dsmcMixedDiffuseSpecularWallRotationPatch
(
    Time& t,
    const polyMesh& mesh,
    dsmcCloud& cloud,
    const dictionary& dict
)
:
    dsmcPatchBoundary(t, mesh, cloud, dict),
    propsDict_(dict.subDict(typeName + "Properties")),
    diffuseFraction_(readScalar(propsDict_.lookup("diffuseFraction")))

{
    writeInTimeDir_ = false;
    writeInCase_ = false;
    measurePropertiesAtWall_ = true;
    
    wallVelocity_ = readScalar(propsDict_.lookup("velocity"));
    wallTemperature_ = readScalar(propsDict_.lookup("temperature"));
    rotationAxis_ = propsDict_.lookup("rotationAxis");
    centrePoint_ = propsDict_.lookup("centrePoint");
    rotationAxis_ /= mag(rotationAxis_);
    
    
    // test

}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

dsmcMixedDiffuseSpecularWallRotationPatch::~dsmcMixedDiffuseSpecularWallRotationPatch()
{}



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
void dsmcMixedDiffuseSpecularWallRotationPatch::initialConfiguration()
{}

void dsmcMixedDiffuseSpecularWallRotationPatch::calculateProperties()
{}

void dsmcMixedDiffuseSpecularWallRotationPatch::controlParticle(dsmcParcel& p, dsmcParcel::trackingData& td)
{
    measurePropertiesBeforeControl(p);

    vector& U = p.U();

    scalar& ERot = p.ERot();
    
    labelList& vibLevel = p.vibLevel();

    label typeId = p.typeId();

    vector nw = p.normal();
    nw /= mag(nw);

    // Normal velocity magnitude
    scalar U_dot_nw = U & nw;

    Random& rndGen(cloud_.rndGen());

    if (diffuseFraction_ > rndGen.scalar01())
    {
        // Diffuse reflection

        // Wall tangential velocity (flow direction)
        vector Ut = U - U_dot_nw*nw;

        while (mag(Ut) < SMALL)
        {
            // If the incident velocity is parallel to the face normal, no
            // tangential direction can be chosen.  Add a perturbation to the
            // incoming velocity and recalculate.

            U = vector
            (
                U.x()*(0.8 + 0.2*rndGen.scalar01()),
                U.y()*(0.8 + 0.2*rndGen.scalar01()),
                U.z()*(0.8 + 0.2*rndGen.scalar01())
            );

            U_dot_nw = U & nw;

            Ut = U - U_dot_nw*nw;
        }

        // Wall tangential unit vector
        vector tw1 = Ut/mag(Ut);

        // Other tangential unit vector
        vector tw2 = nw^tw1;
        
        scalar T = wallTemperature_;

        scalar mass = cloud_.constProps(typeId).mass();

        scalar rotationalDof = cloud_.constProps(typeId).rotationalDegreesOfFreedom();

        scalar vibrationalDof = cloud_.constProps(typeId).vibrationalDegreesOfFreedom();

        U =
            sqrt(physicoChemical::k.value()*T/mass)
           *(
                rndGen.GaussNormal()*tw1
              + rndGen.GaussNormal()*tw2
              - sqrt(-2.0*log(max(1 - rndGen.scalar01(), VSMALL)))*nw
            );
        
        // add wall velocity
        vector uNew = (
                        (p.position() - centrePoint_)/mag((p.position() - centrePoint_))
                      )
                      ^ rotationAxis_;
                      
        uNew /= mag(uNew);
        
        uNew *= wallVelocity_;
        
        U += uNew;
        
        ERot = cloud_.equipartitionRotationalEnergy(T, rotationalDof);

        vibLevel = cloud_.equipartitionVibrationalEnergyLevel(T, vibrationalDof, typeId);
    }
    else
    {
        // Specular reflection

        if (U_dot_nw > 0.0)
        {
            U -= 2.0*U_dot_nw*nw;
        }
        
    }

    measurePropertiesAfterControl(p, 0.0);
}

void dsmcMixedDiffuseSpecularWallRotationPatch::output
(
    const fileName& fixedPathName,
    const fileName& timePath
)
{

}

void dsmcMixedDiffuseSpecularWallRotationPatch::updateProperties(const dictionary& newDict)
{
    //- the main properties should be updated first
    updateBoundaryProperties(newDict);
}



} // End namespace Foam

// ************************************************************************* //
