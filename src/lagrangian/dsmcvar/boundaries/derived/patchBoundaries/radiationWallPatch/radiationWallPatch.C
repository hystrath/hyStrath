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

The wall temperature is a function of the heat flux from incident parcels.

The radiation equation Q = epsilon*sigma*T^4 is used to calculate the wall temperature.

This class gives a single value of temperature for the whole patch.

\*---------------------------------------------------------------------------*/

#include "radiationWallPatch.H"
#include "addToRunTimeSelectionTable.H"
#include "fvc.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(radiationWallPatch, 0);

addToRunTimeSelectionTable(dsmcPatchBoundary, radiationWallPatch, dictionary);



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
radiationWallPatch::radiationWallPatch
(
    Time& t,
    const polyMesh& mesh,
    dsmcCloud& cloud,
    const dictionary& dict
)
:
    dsmcPatchBoundary(t, mesh, cloud, dict),
    propsDict_(dict.subDict(typeName + "Properties")),
    EcTot_(0.0),
    EcTotPrev_(0.0),
    TwallRad_(0.0),
    TwallRadCumul_(0.0),
    epsilonSigma_(0.0),
    thermalCapacity_(0.0),
    alpha_(0.0),
    stepCounter_(0.0)
{
    writeInTimeDir_ = false;
    writeInCase_ = false;
    measurePropertiesAtWall_ = true;

    setProperties();

    TwallRad_ = temperature_;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

radiationWallPatch::~radiationWallPatch()
{}



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
void radiationWallPatch::initialConfiguration()
{}


void radiationWallPatch::calculateProperties()
{
    stepCounter_++;

    const scalar deltaT = mesh_.time().deltaTValue();

    scalar TwallRad = temperature_;

    TwallRad = 
    pow
    (
        alpha_*cloud_.nParticle()*EcTotPrev_/(deltaT*epsilonSigma_*totalPatchSurfaceArea_), 0.25
    ); 

    scalar EcTot = EcTot_;

    // parallel processing of EcTot
    if (Pstream::parRun())
    {
        reduce(EcTot, sumOp<scalar>());
    }

    TwallRadCumul_ += TwallRad 
                      + (1/thermalCapacity_)*deltaT*((cloud_.nParticle()*EcTot/deltaT) 
                      - (totalPatchSurfaceArea_*epsilonSigma_*(pow(TwallRad, 4.0))));


    TwallRad_ = TwallRadCumul_/stepCounter_;

    EcTotPrev_ = EcTot; //stores the EcTot_ value to be used in the next time-step 

    EcTot_ = 0.0;
    
    if(time_.time().outputTime())
    {
        Info << "Temperature at radiation wall patch ( " << patchName_ 
             << "): " << TwallRad_ << endl;
             
        //- reset
        if(resetFieldsAtOutput_)
        {
            stepCounter_ = 0.0;
            TwallRadCumul_ = 0.0;
        }
    }
}

void radiationWallPatch::controlParticle(dsmcParcel& p, dsmcParcel::trackingData& td)
{
    measurePropertiesBeforeControl(p);

    vector& U = p.U();

    scalar& ERot = p.ERot();
    
    labelList& vibLevel = p.vibLevel();

    label typeId = p.typeId();

    // new stuff

    scalar Etrans = 0.5*magSqr(U)*cloud_.constProps(typeId).mass();

    scalar EcTot = Etrans + ERot;

    forAll(vibLevel, i)
    {
        EcTot += vibLevel[i]*physicoChemical::k.value()*cloud_.constProps(typeId).thetaV()[i];
    }
    
    EcTot_ += EcTot;

    //

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

    const scalar& T = TwallRad_;

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

    measurePropertiesAfterControl(p, 0.0);
}

void radiationWallPatch::output
(
    const fileName& fixedPathName,
    const fileName& timePath
)
{

}

void radiationWallPatch::updateProperties(const dictionary& newDict)
{
    //- the main properties should be updated first
    updateBoundaryProperties(newDict);

    propsDict_ = newDict.subDict(typeName + "Properties");

    setProperties();
}

void radiationWallPatch::setProperties()
{
    velocity_ = propsDict_.lookup("velocity");
    temperature_ = readScalar(propsDict_.lookup("temperature"));
    scalar epsilon = readScalar(propsDict_.lookup("emissivity"));
    resetFieldsAtOutput_ = propsDict_.lookup("resetAtOutput");

    epsilonSigma_ = epsilon*physicoChemical::sigma.value();

    thermalCapacity_ = readScalar(propsDict_.lookup("thermalCapacity"));

    alpha_ = readScalar(propsDict_.lookup("energyAccommodationCoeff"));

    if(alpha_ < 0 || alpha_ > 1)
    {

        FatalErrorIn("dsmcRadiationWallFieldPatch::setProperties()")
            << "The value of energyAccommodationCoeff should be between 0 and 1: " 
            << alpha_ << nl 
            << exit(FatalError);
    }
}

} // End namespace Foam

// ************************************************************************* //
