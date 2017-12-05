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

#include "dsmcRadiationWallFieldPatch.H"
#include "addToRunTimeSelectionTable.H"
#include "fvc.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(dsmcRadiationWallFieldPatch, 0);

addToRunTimeSelectionTable(dsmcPatchBoundary, dsmcRadiationWallFieldPatch, dictionary);



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
dsmcRadiationWallFieldPatch::dsmcRadiationWallFieldPatch
(
    Time& t,
    const polyMesh& mesh,
    dsmcCloud& cloud,
    const dictionary& dict
)
:
    dsmcPatchBoundary(t, mesh, cloud, dict),
    propsDict_(dict.subDict(typeName + "Properties")),
    timeDict_(propsDict_.subDict("timeProperties")),
    time_(t, timeDict_),
    time2_(0),
    EcTot_(0.0),
    EcTotPrev_(0.0),
    TwallRad_(0.0),
    TwallRadCumul_(0.0),
    epsilonSigma_(0.0),
    thermalCapacity_(0.0),
    alpha_(0.0)
{
    writeInTimeDir_ = false;
    writeInCase_ = false;
    measurePropertiesAtWall_ = true;

    setProperties();

    TwallRad_ = temperature_;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

dsmcRadiationWallFieldPatch::~dsmcRadiationWallFieldPatch()
{}



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
void dsmcRadiationWallFieldPatch::initialConfiguration()
{}

void dsmcRadiationWallFieldPatch::calculateProperties()
{
    time_++;
    time2_++;

    const scalar deltaT = mesh_.time().deltaTValue();
    
    Info << "Outside loop" << endl;

    if(time2_ > 1)
    {
	  Info << "In loop before TwallRad_" << endl;
    	 TwallRad_ = pow(alpha_*cloud_.nParticle()*EcTotPrev_/(deltaT*epsilonSigma_*patchSurfaceArea_), 0.25); // TwallRad_ here only calculated after first time-step
	 Info << "In loop after TwallRad_" << endl;
    }
    
    Info << "Before TwallRadCumul_" << endl;

    TwallRadCumul_ += TwallRad_ +  (1/thermalCapacity_)*deltaT*((cloud_.nParticle()*EcTot_/deltaT) - (patchSurfaceArea_*epsilonSigma_*(pow(TwallRad_, 4.0))));
    
    Info << "After TwallRadCumul_" << endl;

    if(time_.averagingTime())
    {
        scalar nAvTimeSteps = time_.nAveragingTimeSteps();

        TwallRad_ = TwallRadCumul_/nAvTimeSteps;

        Info << "Temperature at radiation wall patch ( " << patchName_ 
             << "): " << TwallRad_ << endl;

        //- reset 
        if(time_.resetFieldsAtOutput())
        {
            TwallRadCumul_ = 0.0;
        }
	
	EcTotPrev_ = EcTot_; //stores the EcTot_ value to be used in the next time-step 

        EcTot_ = 0.0;
    }


}

void dsmcRadiationWallFieldPatch::controlParticle(dsmcParcel& p, dsmcParcel::trackData& td)
{
    measurePropertiesBeforeControl(p);
    
    Info << "After measurePropertiesBeforeControl" << endl;

    vector& U = p.U();

    scalar& ERot = p.ERot();
    
    scalar& EVib = p.EVib();

    label typeId = p.typeId();


    // new stuff

    scalar Etrans = 0.5*magSqr(U)*cloud_.constProps(typeId).mass();

    scalar EcTot = Etrans + ERot + EVib;
    
    Info << "EcTot" << endl;

   
    EcTot_ += EcTot;

    Info << "EcTot_" << endl;

    vector nw = p.normal();
    nw /= mag(nw);

    // Normal velocity magnitude
    scalar U_dot_nw = U & nw;

    // Wall tangential velocity (flow direction)
    vector Ut = U - U_dot_nw*nw;

    Random& rndGen(cloud_.rndGen());
    
    Info << "Before while loop" << endl;

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
    
    Info << "After while loop" << endl;

    // Wall tangential unit vector
    vector tw1 = Ut/mag(Ut);

    // Other tangential unit vector
    vector tw2 = nw^tw1;
    
    Info << "Before scalar& T" << endl;

    const scalar& T = TwallRad_;
    
    Info << "After scalar& T" << endl;

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
    
    Info << "After velocity generation" << endl;

    ERot = cloud_.equipartitionRotationalEnergy(T, rotationalDof);
    
    EVib = cloud_.equipartitionVibrationalEnergy(T, vibrationalDof, typeId);
    
    Info << "Before measurePropertiesAfterControl" << endl;

    measurePropertiesAfterControl(p);
    
    Info << "After measurePropertiesAfterControl" << endl;
}

void dsmcRadiationWallFieldPatch::output
(
    const fileName& fixedPathName,
    const fileName& timePath
)
{

}

void dsmcRadiationWallFieldPatch::updateProperties(const dictionary& newDict)
{
//     Info << "old dictionary: " << propsDict_  << endl;

    //- the main properties should be updated first
    updateBoundaryProperties(newDict);

    propsDict_ = newDict.subDict(typeName + "Properties");

    setProperties();

}

void dsmcRadiationWallFieldPatch::setProperties()
{
    velocity_ = propsDict_.lookup("velocity");
    temperature_ = readScalar(propsDict_.lookup("temperature"));
    scalar epsilon = readScalar(propsDict_.lookup("emissivity"));
    scalar sigma = readScalar(propsDict_.lookup("stefanBoltzmann"));

    epsilonSigma_ = epsilon* sigma;

    thermalCapacity_ = readScalar(propsDict_.lookup("thermalCapacity"));

    alpha_ = readScalar(propsDict_.lookup("alpha"));

    if(alpha_ < 0 || alpha_ > 1)
    {

        FatalErrorIn("dsmcRadiationWallFieldPatch::setProperties()")
            << "The value of alpha should be between 0 and 1: " 
            << alpha_ << nl 
            << exit(FatalError);
    }

    timeDict_ = propsDict_.subDict("timeProperties");

    if (timeDict_.found("resetAtOutput"))
    {
        time_.resetFieldsAtOutput() = Switch(timeDict_.lookup("resetAtOutput"));
    }
}

} // End namespace Foam

// ************************************************************************* //
