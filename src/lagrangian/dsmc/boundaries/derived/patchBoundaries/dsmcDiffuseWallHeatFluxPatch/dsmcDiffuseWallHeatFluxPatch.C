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

#include "dsmcDiffuseWallHeatFluxPatch.H"
#include "addToRunTimeSelectionTable.H"
#include "fvc.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(dsmcDiffuseWallHeatFluxPatch, 0);

addToRunTimeSelectionTable(dsmcPatchBoundary, dsmcDiffuseWallHeatFluxPatch, dictionary);



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
dsmcDiffuseWallHeatFluxPatch::dsmcDiffuseWallHeatFluxPatch
(
    Time& t,
    const polyMesh& mesh,
    dsmcCloud& cloud,
    const dictionary& dict
)
:
    dsmcPatchBoundary(t, mesh, cloud, dict),
    propsDict_(dict.subDict(typeName + "Properties")),
    faceAreas_(nFaces_,0.0),
    preIE_(0.0),
    postIE_(0.0),
    deltaQ_(nFaces_,0.0),
    mUnUp_(nFaces_,0.0),
    mUn_(nFaces_,0.0),
    uSlip_(nFaces_,0.0),
    resetAtOutputForQ_(true),
    resetCounter_(0.0),
    writeInterval_(0.0),
    firstWrite_(true)
{
    writeInTimeDir_ = true;
    writeInCase_ = true;
    measurePropertiesAtWall_ = true;

    setProperties();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

dsmcDiffuseWallHeatFluxPatch::~dsmcDiffuseWallHeatFluxPatch()
{}



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
void dsmcDiffuseWallHeatFluxPatch::initialConfiguration()
{
    forAll(faceAreas_, f)
    {
        const label& faceI = faces_[f];
        const vector& sF = mesh_.faceAreas()[faceI];
        faceAreas_[f] = mag(sF);
    }
}

void dsmcDiffuseWallHeatFluxPatch::calculateProperties()
{
}

void dsmcDiffuseWallHeatFluxPatch::controlParticle(dsmcParcel& p, dsmcParcel::trackingData& td)
{
    label wppIndex = p.patch(p.face());

    const polyPatch& patch = mesh_.boundaryMesh()[wppIndex];

    label wppLocalFace = patch.whichFace(p.face());

    const dsmcParcel::constantProperties& constProps(cloud_.constProps(p.typeId()));

    const scalar deltaT = mesh_.time().deltaTValue();

    scalar m = constProps.mass();

    preIE_ = 0.5*m*(p.U() & p.U()) + p.ERot();

    forAll(p.vibLevel(), i)
    {
        preIE_ += p.vibLevel()[i]*constProps.thetaV()[i]*physicoChemical::k.value();
    }
    
    measurePropertiesBeforeControl(p);

    vector& U = p.U();

    scalar& ERot = p.ERot();
    
    labelList& vibLevel = p.vibLevel();

    label typeId = p.typeId();

    vector nw = p.normal();
    nw /= mag(nw);

    // Normal velocity magnitude
    scalar U_dot_nw = U & nw;

    // Wall tangential velocity (flow direction)
    vector Ut = U - U_dot_nw*nw;

    Random& rndGen(cloud_.rndGen());

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

    const scalar& T = temperature_;

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

    U += velocity_;

    ERot = cloud_.equipartitionRotationalEnergy(T, rotationalDof);
    
    vibLevel = cloud_.equipartitionVibrationalEnergyLevel(T, vibrationalDof, typeId);

    measurePropertiesAfterControl(p,0.0);

    postIE_ = 0.5*m*(p.U() & p.U()) + p.ERot();

    forAll(p.vibLevel(), i)
    {
        postIE_ += p.vibLevel()[i]*constProps.thetaV()[i]*physicoChemical::k.value();
    }
    
    deltaQ_[wppLocalFace] += cloud_.nParticle()*(preIE_ - postIE_)/(deltaT*faceAreas_[wppLocalFace]);
    
    mUnUp_[wppLocalFace] += ((cloud_.constProps(typeId).mass()/U_dot_nw)*Ut.x());
    
    mUn_[wppLocalFace] += (cloud_.constProps(typeId).mass()/U_dot_nw);
    
    uSlip_[wppLocalFace] = mUnUp_[wppLocalFace]/mUn_[wppLocalFace];

}

void dsmcDiffuseWallHeatFluxPatch::output
(
    const fileName& fixedPathName,
    const fileName& timePath
)
{
    const scalar deltaT = mesh_.time().deltaTValue();
    scalar currentTime = mesh_.time().timeOutputValue();
    
    resetCounter_++;

    if(firstWrite_)
    {
        writeInterval_ = currentTime/(deltaT*resetCounter_);
        firstWrite_ = false;
    }

    if(Pstream::parRun())
    {
        writeTimeData
        (
            timePath,
            patchName()+"_qMean",
            deltaQ_/(writeInterval_*resetCounter_)
        );
        
        writeTimeData
        (
            timePath,
            patchName()+"_uSlip",
            uSlip_
        );
    }
    else
    {
        writeTimeData
        (
            timePath,
            patchName()+"_qMean",
            deltaQ_/(writeInterval_*resetCounter_)
        );	
        
        writeTimeData
        (
            timePath,
            patchName()+"_uSlip",
            uSlip_
        );  
    }

if(resetAtOutputForQ_)
    {
        deltaQ_ = 0.0;
        resetCounter_ = 0.0;
        mUn_ = 0.0;
        mUnUp_ = 0.0;
    }
}

void dsmcDiffuseWallHeatFluxPatch::updateProperties(const dictionary& newDict)
{
    //- the main properties should be updated first
    updateBoundaryProperties(newDict);

    propsDict_ = newDict.subDict(typeName + "Properties");

    setProperties();

}

void dsmcDiffuseWallHeatFluxPatch::setProperties()
{
    velocity_ = propsDict_.lookup("velocity");
    temperature_ = readScalar(propsDict_.lookup("temperature"));
    resetAtOutputForQ_ = propsDict_.lookup("resetAtOutputForQ");
}

} // End namespace Foam

// ************************************************************************* //
