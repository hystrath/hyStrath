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

Cercignani-Lampis-Lord surface interaction model - currently deals with the translation and rotational modes. No vibrational
energy has been dealt with yet and it does not include diffuse scattering with incomplete energy accommodation.

\*---------------------------------------------------------------------------*/

#include "dsmcCLLWallPatch.H"
#include "addToRunTimeSelectionTable.H"
#include "fvc.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(dsmcCLLWallPatch, 0);

addToRunTimeSelectionTable(dsmcPatchBoundary, dsmcCLLWallPatch, dictionary);



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
dsmcCLLWallPatch::dsmcCLLWallPatch
(
    Time& t,
    const polyMesh& mesh,
    dsmcCloud& cloud,
    const dictionary& dict
)
:
    dsmcPatchBoundary(t, mesh, cloud, dict),
    propsDict_(dict.subDict(typeName + "Properties")),
    normalAccommodationCoefficient_(readScalar(propsDict_.lookup("normalAccommodationCoefficient"))),
    tangentialAccommodationCoefficient_(readScalar(propsDict_.lookup("tangentialAccommodationCoefficient"))),
	rotationalEnergyAccommodationCoefficient_(readScalar(propsDict_.lookup("rotationalEnergyAccommodationCoefficient"))),
	vibrationalEnergyAccommodationCoefficient_(readScalar(propsDict_.lookup("vibrationalEnergyAccommodationCoefficient")))
{
    writeInTimeDir_ = false;
    writeInCase_ = false;
    measurePropertiesAtWall_ = true;

    setProperties();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

dsmcCLLWallPatch::~dsmcCLLWallPatch()
{}



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
void dsmcCLLWallPatch::initialConfiguration()
{
    if((normalAccommodationCoefficient_ < VSMALL) && (tangentialAccommodationCoefficient_ < VSMALL))
    {
	measurePropertiesAtWall_ = false;
	//reduces to a specular wall, so no need to measure properties
    }
}

void dsmcCLLWallPatch::calculateProperties()
{

}

void dsmcCLLWallPatch::controlParticle(dsmcParcel& p, dsmcParcel::trackingData& td)
{    
    measurePropertiesBeforeControl(p);

    vector& U = p.U();
    
    label typeId = p.typeId();

    scalar& ERot = p.ERot();
    
//     labelList& vibLevel = p.vibLevel();

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
        
//     scalar vibrationalDof = cloud_.constProps(typeId).vibrationalDegreesOfFreedom();
    
//     scalar characteristicVibrationalTemperature = cloud_.constProps(typeId).thetaV();
    
    const scalar& alphaT = tangentialAccommodationCoefficient_*(2.0 - tangentialAccommodationCoefficient_);
    
    const scalar& alphaN = normalAccommodationCoefficient_;
    
    const scalar& alphaR = rotationalEnergyAccommodationCoefficient_;
    
//     const scalar& alphaV = vibrationalEnergyAccommodationCoefficient_;
    
    scalar mostProbableVelocity = sqrt(2.0*physicoChemical::k.value()*T/mass);
    
    //normalising the incident velocities
    
    vector normalisedTangentialVelocity = Ut/mostProbableVelocity;
    
    scalar normalisedNormalVelocity = U_dot_nw/mostProbableVelocity;
    
    //normal random number components
    
    scalar thetaNormal = 2.0*pi*rndGen.scalar01();
    
    scalar rNormal = sqrt(-alphaN*log(rndGen.scalar01()));

    //tangential random number components

    scalar thetaTangential1 = 2.0*pi*rndGen.scalar01();
    
    scalar rTangential1 = sqrt(-alphaT*log(rndGen.scalar01()));
    
//     scalar thetaTangential2 = 2.0*pi*rndGen.scalar01();
    
//     scalar rTangential2 = sqrt(-alphaT*log(rndGen.scalar01()));

    //selecting the reflected thermal velocities
    
    scalar normalisedIncidentTangentialVelocity1 = mag(normalisedTangentialVelocity);

    scalar um = sqrt(1.0-alphaN)*normalisedNormalVelocity;
    
    scalar normalVelocity = sqrt(
                                    (rNormal*rNormal) 
                                    + (um*um) 
                                    + 2.0*rNormal*um*cos(thetaNormal)
                            );
    
    scalar tangentialVelocity1 = (sqrt(1.0 - alphaT)*mag(normalisedIncidentTangentialVelocity1)
                                + rTangential1*cos(thetaTangential1));
    
    scalar tangentialVelocity2 = rTangential1*sin(thetaTangential1);

    U = 
        mostProbableVelocity
        *(
                tangentialVelocity1*tw1
                + tangentialVelocity2*tw2
                - normalVelocity*nw
        );
            
    if( (p.position().x() > 0.002) && ( p.position().x() < 0.0022) )
    {
        if(Pstream::parRun())
        {
            Pout << "Scattering angle, 0 mm = " << atan(U.y()/U.x()) << endl;
        }
        else
        {
            scalar angle=0;
            if((tangentialVelocity1*tw1).x()>0)
            angle = acos( ( (normalVelocity*nw + tangentialVelocity1*tw1) & tangentialVelocity1*tw1 )/mag((normalVelocity*nw+tangentialVelocity1*tw1) * tangentialVelocity1*tw1) );
            if((tangentialVelocity1*tw1).x()<0)
            angle = acos( ( (normalVelocity*nw - tangentialVelocity1*tw1) & tangentialVelocity1*tw1 )/mag((normalVelocity*nw+tangentialVelocity1*tw1) * tangentialVelocity1*tw1) );
                
            Info << "Scattering angle, 0 mm = " << angle << endl;
        }
    }
    
    if( (p.position().x() > 0.0068) && ( p.position().x() < 0.0072 ) )
    {
        if(Pstream::parRun())
        {
            Pout << "Scattering angle, 5 mm = " << atan(U.y()/U.x()) << endl;
        }
        else
        {
            scalar angle=0;
            if((tangentialVelocity1*tw1).x()>0)
            angle = acos( ( (normalVelocity*nw + tangentialVelocity1*tw1) & tangentialVelocity1*tw1 )/mag((normalVelocity*nw+tangentialVelocity1*tw1) * tangentialVelocity1*tw1) );
            if((tangentialVelocity1*tw1).x()<0)
            angle = acos( ( (normalVelocity*nw - tangentialVelocity1*tw1) & tangentialVelocity1*tw1 )/mag((normalVelocity*nw+tangentialVelocity1*tw1) * tangentialVelocity1*tw1) );
                
            Info << "Scattering angle, 5 mm = " << angle << endl;
        }
    }

    vector uWallNormal = (velocity_ & nw) * nw;
    vector uWallTangential1 = (velocity_ & tw1) * tw1; 
    vector uWallTangential2 = (velocity_ & tw2) * tw2;
    vector UNormal = ((U & nw) * nw) + uWallNormal*alphaN;  
    vector UTangential1 = (U & tw1) * tw1 + uWallTangential1*alphaT;
    vector UTangential2 = (U & tw2) * tw2 + uWallTangential2*alphaT;
    
    U = UNormal + UTangential1 + UTangential2;
    
    //selecting rotational energy, this is Lord's extension to rotational degrees of freedom
    
    if(rotationalDof == 2)
    {
        scalar om = sqrt( (ERot*(1.0 - alphaR)) / (physicoChemical::k.value()*T));
        scalar rRot = sqrt(-alphaR*(log(max(1.0 - rndGen.scalar01(), VSMALL))));
        scalar thetaRot = 2.0*pi*rndGen.scalar01();
        ERot = physicoChemical::k.value()*T*((rRot*rRot) + (om*om) + (2.0*rRot*om*cos(thetaRot)));
    }
    
    if(rotationalDof == 3) //polyatomic case, see Bird's DSMC2.FOR code
    {
        scalar X = 0.0; 
        scalar A = 0.0;
        
        do
        {
                X = 4.0*rndGen.scalar01();
                A = 2.7182818*X*X*exp(-(X*X));
        } while (A < rndGen.scalar01());
        
        scalar om = sqrt( (ERot*(1.0 - alphaR)) / (physicoChemical::k.value()*T));
        scalar rRot = sqrt(-alphaR)*X;
        scalar thetaRot = 2.0*rndGen.scalar01() - 1.0;
        ERot = physicoChemical::k.value()*T*((rRot*rRot) + (om*om) + (2.0*rRot*om*cos(thetaRot)));
    }
	
// 	//selecting vibrational energy, this is lord's extension to vibrational degrees of freedom
// 	
// 	if(vibrationalDof > VSMALL)
// 	{
// 		scalar EVibStar = -log(1.0 - rndGen.scalar01() + rndGen.scalar01()*exp(-characteristicVibrationalTemperature*physicoChemical::k.value())); 
// 		
// 		EVib += EVibStar;
// 		
// 		scalar vm = sqrt(1.0-alphaV)*sqrt(EVib);
// 		
// 		scalar thetaVib = 2.0*pi*rndGen.scalar01();
// 		
// 		scalar rVib = sqrt(-alphaV*log(rndGen.scalar01()));
// 		
// 		EVib = (
// 					(rVib*rVib) 
// 					+ (vm*vm) 
// 					+ 2.0*rVib*vm*cos(thetaVib)
// 				);
// 		
// 		label iPost = EVib / (characteristicVibrationalTemperature*physicoChemical::k.value()); //post interaction vibrational quantum level
// 		
// 		EVib = iPost * characteristicVibrationalTemperature*physicoChemical::k.value();
// 	}

// 	EVib = cloud_.equipartitionVibrationalEnergy(T, vibrationalDof, typeId);


    measurePropertiesAfterControl(p, 0.0);
}

void dsmcCLLWallPatch::output
(
    const fileName& fixedPathName,
    const fileName& timePath
)
{

}

void dsmcCLLWallPatch::updateProperties(const dictionary& newDict)
{
//     Info << "old dictionary: " << propsDict_  << endl;

    //- the main properties should be updated first
    updateBoundaryProperties(newDict);

    propsDict_ = newDict.subDict(typeName + "Properties");

    setProperties();

}

void dsmcCLLWallPatch::setProperties()
{
    velocity_ = propsDict_.lookup("velocity");
    temperature_ = readScalar(propsDict_.lookup("temperature"));
}

} // End namespace Foam

// ************************************************************************* //
