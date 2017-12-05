/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2009-2010 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

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

\*---------------------------------------------------------------------------*/

#include "LarsenBorgnakkeVariableSoftSphere.H"
#include "constants.H"
#include "addToRunTimeSelectionTable.H"

using namespace Foam::constant::mathematical;

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(LarsenBorgnakkeVariableSoftSphere, 0);
    addToRunTimeSelectionTable(BinaryCollisionModel, 
                               LarsenBorgnakkeVariableSoftSphere, dictionary);
};

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


Foam::LarsenBorgnakkeVariableSoftSphere::LarsenBorgnakkeVariableSoftSphere
(
    const dictionary& dict,
    dsmcCloud& cloud
)
:
    BinaryCollisionModel(dict, cloud),
    coeffDict_(dict.subDict(typeName + "Coeffs")),
    Tref_(readScalar(coeffDict_.lookup("Tref"))),
    rotationalRelaxationCollisionNumber_
    (
        readScalar(coeffDict_.lookup("rotationalRelaxationCollisionNumber"))
    ),
    electronicRelaxationCollisionNumber_
    (
        readScalar(coeffDict_.lookup("electronicRelaxationCollisionNumber"))
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //


Foam::LarsenBorgnakkeVariableSoftSphere::~LarsenBorgnakkeVariableSoftSphere()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


bool Foam::LarsenBorgnakkeVariableSoftSphere::active() const
{
    return true;
}



Foam::scalar Foam::LarsenBorgnakkeVariableSoftSphere::sigmaTcR
(
    const dsmcParcel& pP,
    const dsmcParcel& pQ
) const
{
//     const CloudType& cloud(this->owner());

    label typeIdP = pP.typeId();
    label typeIdQ = pQ.typeId();

    scalar dPQ =
        0.5
       *(
            cloud_.constProps(typeIdP).d()
          + cloud_.constProps(typeIdQ).d()
        );

    scalar omegaPQ =
        0.5
       *(
            cloud_.constProps(typeIdP).omega()
          + cloud_.constProps(typeIdQ).omega()
        );

    scalar cR = mag(pP.U() - pQ.U());

    if (cR < VSMALL)
    {
        return 0;
    }

    scalar mP = cloud_.constProps(typeIdP).mass();

    scalar mQ = cloud_.constProps(typeIdQ).mass();

    scalar mR = mP*mQ/(mP + mQ);

    // calculating cross section = pi*dPQ^2, where dPQ is from Bird, eq. 4.79
    scalar sigmaTPQ =
        pi*dPQ*dPQ
       *pow(2.0*physicoChemical::k.value()*Tref_/(mR*cR*cR), omegaPQ - 0.5)
       /exp(Foam::lgamma(2.5 - omegaPQ));

    return sigmaTPQ*cR;
}



void Foam::LarsenBorgnakkeVariableSoftSphere::collide
(
    dsmcParcel& pP,
    dsmcParcel& pQ,
    label& cellI
)
{

    label typeIdP = pP.typeId();
    label typeIdQ = pQ.typeId();
    vector& UP = pP.U();
    vector& UQ = pQ.U();
    scalar& ERotP = pP.ERot();
    scalar& ERotQ = pQ.ERot();
    label& ELevelP = pP.ELevel();
    label& ELevelQ = pQ.ELevel();
    labelList& vibLevelP = pP.vibLevel();
    labelList& vibLevelQ = pQ.vibLevel();
    
    scalar alphaPQ = 
    0.5*(
        cloud_.constProps(typeIdP).alpha()
        + cloud_.constProps(typeIdQ).alpha()
    );
    
    
    scalar collisionSeparation = sqrt(
            sqr(pP.position().x() - pQ.position().x()) +
            sqr(pP.position().y() - pQ.position().y())
    );
    
    cloud_.cellPropMeasurements().collisionSeparation()[cellI] += 
                                                        collisionSeparation;
    cloud_.cellPropMeasurements().nColls()[cellI]++;

    Random& rndGen(cloud_.rndGen());
    
    //Larsen Borgnakke rotational energy redistribution part.  Using the serial
    // application of the LB method, as per the INELRS subroutine in Bird's
    //DSMC0R.FOR

    scalar preCollisionERotP = ERotP;
    scalar preCollisionERotQ = ERotQ;
    
    scalarList preCollisionEVibP(vibLevelP.size(),0.0);
    scalarList preCollisionEVibQ(vibLevelQ.size(),0.0);
    
    scalar vibrationalDofP = 
                cloud_.constProps(typeIdP).vibrationalDegreesOfFreedom();
    scalar vibrationalDofQ = 
                cloud_.constProps(typeIdQ).vibrationalDegreesOfFreedom();
                
    if(vibrationalDofP > 0)
    {
        forAll(vibLevelP, i)
        {
            preCollisionEVibP[i] =  vibLevelP[i]*
                                    cloud_.constProps(typeIdP).thetaV()[i]
                                    *physicoChemical::k.value();
        }
    }

    if(vibrationalDofQ > 0)
    {
        forAll(vibLevelQ, i)
        {
            preCollisionEVibQ[i] =  vibLevelQ[i]*
                                    cloud_.constProps(typeIdQ).thetaV()[i]
                                    *physicoChemical::k.value();
        }
    }

    scalar preCollisionEEleP = 
                cloud_.constProps(typeIdP).electronicEnergyList()[ELevelP];
    scalar preCollisionEEleQ = 
                cloud_.constProps(typeIdQ).electronicEnergyList()[ELevelQ];

    scalar rotationalDofP = 
                    cloud_.constProps(typeIdP).rotationalDegreesOfFreedom();

    scalar rotationalDofQ = 
                    cloud_.constProps(typeIdQ).rotationalDegreesOfFreedom();
                    
    label jMaxP = cloud_.constProps(typeIdP).numberOfElectronicLevels();    
    label jMaxQ = cloud_.constProps(typeIdQ).numberOfElectronicLevels();
    
    List<scalar> EElistP = cloud_.constProps(typeIdP).electronicEnergyList();   
 
    List<scalar> EElistQ = cloud_.constProps(typeIdQ).electronicEnergyList();
   
    List<label> gListP = cloud_.constProps(typeIdP).degeneracyList();    
    List<label> gListQ = cloud_.constProps(typeIdQ).degeneracyList();  
    
    scalarList thetaVP = cloud_.constProps(typeIdP).thetaV();  
    scalarList thetaVQ = cloud_.constProps(typeIdQ).thetaV();
    
    scalarList thetaDP = cloud_.constProps(typeIdP).thetaD();
    scalarList thetaDQ = cloud_.constProps(typeIdQ).thetaD();
    
    scalarList ZrefP = cloud_.constProps(typeIdP).Zref();
    scalarList ZrefQ = cloud_.constProps(typeIdQ).Zref();
    
    scalarList refTempZvP = cloud_.constProps(typeIdP).TrefZv();
    scalarList refTempZvQ = cloud_.constProps(typeIdQ).TrefZv();

    scalar omegaPQ =
        0.5
       *(
            cloud_.constProps(typeIdP).omega()
          + cloud_.constProps(typeIdQ).omega()
        );

    scalar mP = cloud_.constProps(typeIdP).mass();
    scalar mQ = cloud_.constProps(typeIdQ).mass();

    scalar mR = mP*mQ/(mP + mQ);
    vector Ucm = (mP*UP + mQ*UQ)/(mP + mQ);
    scalar cRsqr = magSqr(UP - UQ);

    scalar translationalEnergy = 0.5*mR*cRsqr;

    scalar ChiB = 2.5 - omegaPQ;
    
    scalar inverseRotationalCollisionNumber = 
                                    1.0/rotationalRelaxationCollisionNumber_;
    scalar inverseElectronicCollisionNumber = 
                                    1.0/electronicRelaxationCollisionNumber_;
                                    
    // Larsen Borgnakke rotational energy redistribution part.  Using the 
    // serial application of the LB method, as per the INELRS subroutine in 
    // Bird's DSMC0R.FOR
    
    if(inverseElectronicCollisionNumber > rndGen.scalar01())
    { 

        // collision energy of particle P = relative translational energy + 
        // pre-collision electronic energy
        
        scalar EcP = translationalEnergy + preCollisionEEleP;
        
        label postCollisionELevel = cloud_.postCollisionElectronicEnergyLevel
                        (
                            EcP,
                            jMaxP,
                            omegaPQ,
                            EElistP,
                            gListP
                        );
                        
        ELevelP = postCollisionELevel;
        
        // relative translational energy after electronic exchange
        translationalEnergy = EcP - EElistP[ELevelP];
    }
            
    if(vibrationalDofP > VSMALL)
    {
        forAll(vibLevelP, i)
        {
            // collision energy of particle P = relative translational energy + 
            // pre-collision vibrational energy
            scalar EcP = translationalEnergy + preCollisionEVibP[i]; 

            // - maximum possible quantum level (equation 3, Bird 2010)
            label iMaxP = (EcP / (physicoChemical::k.value()*thetaVP[i])); 

            if(iMaxP > SMALL)
            {       
                vibLevelP[i] = cloud_.postCollisionVibrationalEnergyLevel
                        (
                            false,
                            vibLevelP[i],
                            iMaxP,
                            thetaVP[i],
                            thetaDP[i],
                            refTempZvP[i],
                            omegaPQ,
                            ZrefP[i],
                            EcP
                        );
                        
                translationalEnergy = EcP 
                                        - (vibLevelP[i]
                                        *cloud_.constProps(typeIdP).thetaV()[i]
                                        *physicoChemical::k.value());
            }
        }
    }

    if (rotationalDofP > 0)
    {
        if (inverseRotationalCollisionNumber > rndGen.scalar01())
        {
            scalar EcP = translationalEnergy + preCollisionERotP;
            
            scalar energyRatio = 
                    cloud_.postCollisionRotationalEnergy(rotationalDofP,ChiB);

            ERotP = energyRatio*EcP;
        
            translationalEnergy = EcP - ERotP;
        }
    }
    
    if(inverseElectronicCollisionNumber > rndGen.scalar01())
    {

        // collision energy of particle Q = relative translational energy + 
        // pre-collision electronic energy
        
        scalar EcQ = translationalEnergy + preCollisionEEleQ; 
        
        label postCollisionELevel = cloud_.postCollisionElectronicEnergyLevel
                        (
                            EcQ,
                            jMaxQ,
                            omegaPQ,
                            EElistQ,
                            gListQ
                        );
                        
        ELevelQ = postCollisionELevel;

        // relative translational energy after electronic exchange
        translationalEnergy = EcQ - EElistQ[ELevelQ];
    }
              
    if(vibrationalDofQ > VSMALL)
    {
        forAll(vibLevelQ, i)
        {
            // collision energy of particle Q = relative translational energy + 
            // pre-collision vibrational energy
            scalar EcQ = translationalEnergy + preCollisionEVibQ[i]; 

            // - maximum possible quantum level (equation 3, Bird 2010)
            label iMaxQ = (EcQ / (physicoChemical::k.value()*thetaVQ[i])); 

            if(iMaxQ > SMALL)
            {       
                vibLevelQ[i] = cloud_.postCollisionVibrationalEnergyLevel
                        (
                            false,
                            vibLevelQ[i],
                            iMaxQ,
                            thetaVQ[i],
                            thetaDQ[i],
                            refTempZvQ[i],
                            omegaPQ,
                            ZrefQ[i],
                            EcQ
                        );
                        
                translationalEnergy = EcQ 
                                        - (vibLevelQ[i]
                                        *cloud_.constProps(typeIdQ).thetaV()[i]
                                        *physicoChemical::k.value());
            }
        }
    }

    if (rotationalDofQ > 0)
    {
        if (inverseRotationalCollisionNumber > rndGen.scalar01())
        {
            scalar EcQ = translationalEnergy + preCollisionERotQ;
            
            scalar energyRatio = 
                    cloud_.postCollisionRotationalEnergy(rotationalDofQ,ChiB);

            ERotQ = energyRatio*EcQ;
        
            translationalEnergy = EcQ - ERotQ;
        }
    }

    // Rescale the translational energy
    scalar A = sqrt(2.0*translationalEnergy/mR);
    
    scalar cR = mag(UP - UQ);
    
    vector cRComponents = (UP - UQ) * (A/cR);
    
    cR = A;
   
    scalar cosTheta = (2.0*pow(rndGen.scalar01(),(1.0/alphaPQ))) - 1.0;

    scalar sinTheta = sqrt(1.0 - cosTheta*cosTheta);

    scalar phi = twoPi*rndGen.scalar01();
    
    scalar D = sqrt(pow(cRComponents.y(),2.0) + pow(cRComponents.z(),2.0));

    vector postCollisionRelU =
        vector
        (
            (cosTheta*cRComponents.x()) + sinTheta*sin(phi)*D,
            (cosTheta*cRComponents.y()) + sinTheta
                *(cR*cRComponents.z()*cos(phi) 
                    - cRComponents.x()*cRComponents.y()
                    *sin(phi))/D,
            (cosTheta*cRComponents.z()) - sinTheta
                *(cR*cRComponents.y()*cos(phi) 
                + cRComponents.x()*cRComponents.z()
                *sin(phi))/D
        );

    UP = Ucm + postCollisionRelU*mQ/(mP + mQ);

    UQ = Ucm - postCollisionRelU*mP/(mP + mQ);
    
    label classificationP = pP.classification();
    label classificationQ = pQ.classification();
    
    //- class I molecule changes to class
    //- III molecule when it collides with either class II or class III
    //- molecules.
    
    if(classificationP == 0 && classificationQ == 1)
    {
        pP.classification() = 2;
    }
    
    if(classificationQ == 0 && classificationP == 1)
    {
        pQ.classification() = 2;
    }
    
    if(classificationP == 0 && classificationQ == 2)
    {
        pP.classification() = 2;
    }
    
    if(classificationQ == 0 && classificationP == 2)
    {
        pQ.classification() = 2;
    }
}


const Foam::dictionary& 
Foam::LarsenBorgnakkeVariableSoftSphere::coeffDict() const
{
    return coeffDict_;
}


// ************************************************************************* //
