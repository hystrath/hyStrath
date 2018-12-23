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

#include "dsmcDiffuseWallZeroGradMultiTPatch.H"
#include "addToRunTimeSelectionTable.H"
#include "fvc.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(dsmcDiffuseWallZeroGradMultiTPatch, 0);

addToRunTimeSelectionTable(dsmcPatchBoundary, dsmcDiffuseWallZeroGradMultiTPatch, dictionary);



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
dsmcDiffuseWallZeroGradMultiTPatch::dsmcDiffuseWallZeroGradMultiTPatch
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

    setProperties();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

dsmcDiffuseWallZeroGradMultiTPatch::~dsmcDiffuseWallZeroGradMultiTPatch()
{}



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
void dsmcDiffuseWallZeroGradMultiTPatch::initialConfiguration()
{
    
}

void dsmcDiffuseWallZeroGradMultiTPatch::calculateProperties()
{

}

void dsmcDiffuseWallZeroGradMultiTPatch::controlParticle(dsmcParcel& p, dsmcParcel::trackingData& td)
{
    measurePropertiesBeforeControl(p);
    
//    scalar currentTime = cloud_.mesh().time().value();
    
//     Info << "currentTime = " << currentTime << endl;

    vector& U = p.U();

    scalar& ERot = p.ERot();
    
    labelList& vibLevel = p.vibLevel();
    
    label& ELevel = p.ELevel();

    const label& typeId = p.typeId();

    vector nw = p.normal();
    nw /= mag(nw);

    // Normal velocity magnitude
    scalar U_dot_nw = U & nw;

    // Wall tangential velocity (flow direction)
    vector Ut = U - U_dot_nw*nw;

    Random& rndGen(cloud_.rndGen());

    const scalar& Ttra = temperatureT_;
    const scalar& Trot = temperatureR_;
    const scalar& Tvib = temperatureV_;
    const scalar& Tel = temperatureE_;

    scalar mass = cloud_.constProps(typeId).mass();

    scalar rotationalDof = cloud_.constProps(typeId).rotationalDegreesOfFreedom();
    
    scalar vibrationalDof = cloud_.constProps(typeId).nVibrationalModes();

    if (zeroGradTt_)
    {
        if (U_dot_nw > 0.0)
        {
            U -= 2.0*U_dot_nw*nw;
        }
    }
    else
    {
        while (mag(Ut) < SMALL)
        {
            // If the incident velocity is parallel to the face normal, no
            // tangential direction can be chosen.  Add a perturbation to the
            // incoming velocity and recalculate.

            U = vector
            (
                U.x()*(0.8 + 0.2*rndGen.sample01<scalar>()),
                U.y()*(0.8 + 0.2*rndGen.sample01<scalar>()),
                U.z()*(0.8 + 0.2*rndGen.sample01<scalar>())
            );

            U_dot_nw = U & nw;

            Ut = U - U_dot_nw*nw;
        }

        // Wall tangential unit vector
        vector tw1 = Ut/mag(Ut);

        // Other tangential unit vector
        vector tw2 = nw^tw1;

        U = sqrt(physicoChemical::k.value()*Ttra/mass) * (
            rndGen.GaussNormal<scalar>() * tw1
          + rndGen.GaussNormal<scalar>() * tw2
          - sqrt(-2.0*log(max(1 - rndGen.sample01<scalar>(), VSMALL))) * nw
            );

        U += velocity_;
    }


    if(not zeroGradTr_)
    {
        ERot = cloud_.equipartitionRotationalEnergy(Trot, rotationalDof);
    }

    if(not zeroGradTv_)
    {
        vibLevel = cloud_.equipartitionVibrationalEnergyLevel
            (
                Tvib, 
                vibrationalDof, 
                typeId
            );
    }
    
    if(not zeroGradTe_)
    {
        ELevel = cloud_.equipartitionElectronicLevel
        (
            Tel,
            cloud_.constProps(typeId).electronicDegeneracyList(),
            cloud_.constProps(typeId).electronicEnergyList(),
            typeId
        );
    }


    measurePropertiesAfterControl(p, 0.0);
}

void dsmcDiffuseWallZeroGradMultiTPatch::output
(
    const fileName& fixedPathName,
    const fileName& timePath
)
{
}


void dsmcDiffuseWallZeroGradMultiTPatch::updateProperties(const dictionary& newDict)
{
    //- the main properties should be updated first
    updateBoundaryProperties(newDict);

    propsDict_ = newDict.subDict(typeName + "Properties");

    setProperties();

}

void dsmcDiffuseWallZeroGradMultiTPatch::setProperties()
{
    velocity_ = propsDict_.lookup("velocity");
    temperatureT_ = readScalar(propsDict_.lookup("translationalT"));
    temperatureR_ = readScalar(propsDict_.lookup("rotationalT"));
    temperatureV_ = readScalar(propsDict_.lookup("vibrationalT"));
    temperatureE_ = readScalar(propsDict_.lookup("electronicT"));
    zeroGradTt_ = propsDict_.lookupOrDefault<bool>("zeroGradTt", false);
    zeroGradTr_ = propsDict_.lookupOrDefault<bool>("zeroGradTr", false);
    zeroGradTv_ = propsDict_.lookupOrDefault<bool>("zeroGradTv", false);
    zeroGradTe_ = propsDict_.lookupOrDefault<bool>("zeroGradTe", false);
    Info << zeroGradTt_ << zeroGradTr_ << zeroGradTv_ << endl;
}

} // End namespace Foam

// ************************************************************************* //

