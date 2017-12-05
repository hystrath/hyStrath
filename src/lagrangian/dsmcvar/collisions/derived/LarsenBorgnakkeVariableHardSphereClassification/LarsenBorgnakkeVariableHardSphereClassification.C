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

#include "LarsenBorgnakkeVariableHardSphereClassification.H"
#include "constants.H"
#include "addToRunTimeSelectionTable.H"

using namespace Foam::constant::mathematical;

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(LarsenBorgnakkeVariableHardSphereClassification, 0);
    addToRunTimeSelectionTable(BinaryCollisionModel, LarsenBorgnakkeVariableHardSphereClassification, dictionary);
};

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

//template <class CloudType>
Foam::LarsenBorgnakkeVariableHardSphereClassification::LarsenBorgnakkeVariableHardSphereClassification
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
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

//template <class CloudType>
Foam::LarsenBorgnakkeVariableHardSphereClassification::~LarsenBorgnakkeVariableHardSphereClassification()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::LarsenBorgnakkeVariableHardSphereClassification::active() const
{
    return true;
}

//template <class CloudType>
Foam::scalar Foam::LarsenBorgnakkeVariableHardSphereClassification::sigmaTcR
(
    const dsmcParcel& pP,
    const dsmcParcel& pQ
) const
{
    //const CloudType& cloud(this->owner());
    
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


//template <class CloudType>
void Foam::LarsenBorgnakkeVariableHardSphereClassification::collide
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
    scalarList EVibP(pP.vibLevel().size(), 0.0);
    scalarList EVibQ(pQ.vibLevel().size(), 0.0);
    label& ELevelP = pP.ELevel();
    label& ELevelQ = pQ.ELevel();
    labelList& vibLevelP = pP.vibLevel();
    labelList& vibLevelQ = pQ.vibLevel();
    
    forAll(EVibP, i)
    {
        EVibP[i] = pP.vibLevel()[i]*cloud_.constProps(typeIdP).thetaV()[i]*physicoChemical::k.value();
    }
    
    forAll(EVibQ, i)
    {
        EVibQ[i] = pQ.vibLevel()[i]*cloud_.constProps(typeIdQ).thetaV()[i]*physicoChemical::k.value();
    }
    
    scalar collisionSeparation = sqrt(
            sqr(pP.position().x() - pQ.position().x()) +
            sqr(pP.position().y() - pQ.position().y())
    );
    
    cloud_.cellPropMeasurements().collisionSeparation()[cellI] += collisionSeparation;
    cloud_.cellPropMeasurements().nColls()[cellI]++;

    Random& rndGen(cloud_.rndGen());
    
    
  //   VIBRATIONAL ENERGY EXCHANGE - QUANTUM-KINETIC MODEL
    
    scalar preCollisionERotP = ERotP;
    scalar preCollisionERotQ = ERotQ;
    
    scalarList preCollisionEVibP = EVibP;
    scalarList preCollisionEVibQ = EVibQ;
    
//     scalar preCollisionEEleP = cloud_.constProps(typeIdP).electronicEnergyList()[ELevelP];
//     scalar preCollisionEEleQ = cloud_.constProps(typeIdQ).electronicEnergyList()[ELevelQ];

    scalar rotationalDofP = cloud_.constProps(typeIdP).rotationalDegreesOfFreedom();
    scalar rotationalDofQ = cloud_.constProps(typeIdQ).rotationalDegreesOfFreedom();
    
    scalar vibrationalDofP = cloud_.constProps(typeIdP).vibrationalDegreesOfFreedom();
    scalar vibrationalDofQ = cloud_.constProps(typeIdQ).vibrationalDegreesOfFreedom();
    
    label jMaxP = cloud_.constProps(typeIdP).numberOfElectronicLevels();    
    label jMaxQ = cloud_.constProps(typeIdQ).numberOfElectronicLevels();
    
//     List<scalar> EElistP = cloud_.constProps(typeIdP).electronicEnergyList();    
//     List<scalar> EElistQ = cloud_.constProps(typeIdQ).electronicEnergyList();
   
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
    
    scalar inverseRotationalCollisionNumber = 1.0/rotationalRelaxationCollisionNumber_;
    
    // Larsen Borgnakke rotational energy redistribution part.  Using the serial
    // application of the LB method, as per the INELRS subroutine in Bird's
    // DSMC0R.FOR
            
    if(vibrationalDofP > VSMALL)
    {
        forAll(EVibP, i)
        {
            // collision energy of particle P = relative translational energy + pre-collision vibrational energy
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
                        
                translationalEnergy = EcP - (vibLevelP[i]*cloud_.constProps(typeIdP).thetaV()[i]*physicoChemical::k.value());
            }
        }
    }
        
    if (rotationalDofP > VSMALL)
    {
//         scalar particleProbabilityP = 
//             ((zeta_T + 2.0*rotationalDofP)/(2.0*rotationalDofP))
//             *(
//                 1.0 - sqrt(
//                             1.0 - (rotationalDofP/zeta_T)
//                             *((zeta_T+rotationalDofP)/(zeta_T+2.0*rotationalDofP))
//                             *(4.0/rotationalRelaxationCollisionNumber_)
//                           )
//              );
            
//         Info << "particleProbabilityP = " << particleProbabilityP << endl;
        
        if (inverseRotationalCollisionNumber /*particleProbabilityP*/ > rndGen.scalar01())
        {
            scalar EcP = translationalEnergy + preCollisionERotP;
            
            scalar energyRatio = cloud_.postCollisionRotationalEnergy(rotationalDofP,ChiB);

            ERotP = energyRatio*EcP;
        
            translationalEnergy = EcP - ERotP;
        }
    }   

      
    if(vibrationalDofQ > VSMALL)
    {
        forAll(EVibQ, i)
        {
            // collision energy of particle Q = relative translational energy + pre-collision vibrational energy
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
                        
                translationalEnergy = EcQ - (vibLevelQ[i]*cloud_.constProps(typeIdQ).thetaV()[i]*physicoChemical::k.value());
            }
        }
    }
        
    if (rotationalDofQ > VSMALL)
    {
//         scalar particleProbabilityQ = 
//             ((zeta_T + 2.0*rotationalDofQ)/(2.0*rotationalDofQ))
//             *(
//                 1.0 - sqrt(
//                             1.0 - (rotationalDofQ/zeta_T)
//                             *((zeta_T+rotationalDofQ)/(zeta_T+2.0*rotationalDofQ))
//                             *(4.0/rotationalRelaxationCollisionNumber_)
//                           )
//              );
            
//         Info << "particleProbabilityQ = " << particleProbabilityQ << endl;
        
        if (inverseRotationalCollisionNumber /*particleProbabilityQ*/ > rndGen.scalar01())
        {
            scalar EcQ = translationalEnergy + preCollisionERotQ;
            
            scalar energyRatio = cloud_.postCollisionRotationalEnergy(rotationalDofQ,ChiB);

            ERotQ = energyRatio*EcQ;
        
            translationalEnergy = EcQ - ERotQ;
        }
    }

    // Rescale the translational energy
    scalar cR = sqrt((2.0*translationalEnergy)/mR);

    // Variable Hard Sphere collision part

    scalar cosTheta = 2.0*rndGen.scalar01() - 1.0;

    scalar sinTheta = sqrt(1.0 - cosTheta*cosTheta);

    scalar phi = 2.0*pi*rndGen.scalar01();

    vector postCollisionRelU =
        cR
       *vector
        (
            cosTheta,
            sinTheta*cos(phi),
            sinTheta*sin(phi)
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

const Foam::dictionary& Foam::LarsenBorgnakkeVariableHardSphereClassification::coeffDict() const
{
    return coeffDict_;
}


// ************************************************************************* //
