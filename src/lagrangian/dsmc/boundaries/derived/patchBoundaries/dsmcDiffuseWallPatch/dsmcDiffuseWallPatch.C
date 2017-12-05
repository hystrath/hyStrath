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

#include "dsmcDiffuseWallPatch.H"
#include "addToRunTimeSelectionTable.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(dsmcDiffuseWallPatch, 0);

addToRunTimeSelectionTable
(
    dsmcPatchBoundary, 
    dsmcDiffuseWallPatch, 
    dictionary
);


// * * * * * * * * * * * Protected member functions  * * * * * * * * * * * * //

void dsmcDiffuseWallPatch::setProperties()
{
    velocity_ = propsDict_.lookup("velocity");
    
    if(propsDict_.found("groundLevelTemperature"))
    {
        temperature_ = readScalar(propsDict_.lookup("groundLevelTemperature"));
    }
    else
    {
        temperature_ = readScalar(propsDict_.lookup("temperature"));
    }
    
    formationLevelTemperature_ = 
        propsDict_.lookupOrDefault<scalar>("formationLevelTemperature", temperature_);
}


void dsmcDiffuseWallPatch::performDiffuseReflection
(
    dsmcParcel& p,
    const scalar& localTemperature,
    const vector& localVelocity
)
{
    scalar T = getLocalTemperature(p.position()[depthAxis_]);
    
    if(localTemperature != 0)
    {
        T = localTemperature;
    }
    
    vector& U = p.U();

    scalar& ERot = p.ERot();
    
    labelList& vibLevel = p.vibLevel();
    
    label& ELevel = p.ELevel();

    const label& typeId = p.typeId();
    
    //- Wall unit normal vector and wall unit tangential vectors
    vector nw, tw1, tw2 = vector::zero;

    dsmcPatchBoundary::calculateWallUnitVectors(p, nw, tw1, tw2);

    const scalar mass = cloud_.constProps(typeId).mass();

    scalar rotationalDof = 
        cloud_.constProps(typeId).rotationalDegreesOfFreedom();
    
    scalar vibrationalDof = 
        cloud_.constProps(typeId).vibrationalDegreesOfFreedom();

    Random& rndGen = cloud_.rndGen();

    U = sqrt(physicoChemical::k.value()*T/mass)
       *(
            rndGen.GaussNormal()*tw1
          + rndGen.GaussNormal()*tw2
          - sqrt(-2.0*log(max(1 - rndGen.scalar01(), VSMALL)))*nw
        );

       
    ERot = cloud_.equipartitionRotationalEnergy(T, rotationalDof);

    
    vibLevel = 
        cloud_.equipartitionVibrationalEnergyLevel(T, vibrationalDof, typeId);
   
    
    ELevel = cloud_.equipartitionElectronicLevel
        (
            T,
            cloud_.constProps(typeId).degeneracyList(),
            cloud_.constProps(typeId).electronicEnergyList(),
            typeId
        );   
    
    if(localVelocity != vector::zero)
    {
        U += localVelocity;
    }
    else
    {
        U += velocity_;
    }
}


scalar dsmcDiffuseWallPatch::getLocalTemperature
(
    const scalar& depthPosition
) const
{
    return temperature_ + (depthPosition - maxDepth_)
        *(temperature_ - formationLevelTemperature_)/lengthPatch_;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
dsmcDiffuseWallPatch::dsmcDiffuseWallPatch
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
    
    depthAxis_ = 1;
    
    word depthAxisChar = propsDict_.lookupOrDefault<word>("depthAxis", "y");
    
    if (depthAxisChar == "x")
    {
        depthAxis_ = 0;
    }
    else if (depthAxisChar == "z")
    {
        depthAxis_ = 2;
    }
    
    // Here "mesh" refers to the polyMesh (patch) so the min and max positions
    // should be that of the patch bounding box
    maxDepth_ = mesh.bounds().max().component(depthAxis_);
    lengthPatch_ = maxDepth_ - mesh.bounds().min().component(depthAxis_);

    setProperties();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

dsmcDiffuseWallPatch::~dsmcDiffuseWallPatch()
{}



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void dsmcDiffuseWallPatch::initialConfiguration()
{}


void dsmcDiffuseWallPatch::calculateProperties()
{}


void dsmcDiffuseWallPatch::controlParticle
(
    dsmcParcel& p, 
    dsmcParcel::trackingData& td
)
{
    measurePropertiesBeforeControl(p);
    
    performDiffuseReflection(p);
    
    measurePropertiesAfterControl(p);
}


void dsmcDiffuseWallPatch::output
(
    const fileName& fixedPathName,
    const fileName& timePath
)
{}


void dsmcDiffuseWallPatch::updateProperties(const dictionary& newDict)
{
    //- the main properties should be updated first
    updateBoundaryProperties(newDict);

    propsDict_ = newDict.subDict(typeName + "Properties");

    setProperties();
}

} // End namespace Foam

// ************************************************************************* //
