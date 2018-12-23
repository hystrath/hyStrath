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

#include "LarsenBorgnakkeVariableHardSphere.H"
#include "constants.H"
#include "addToRunTimeSelectionTable.H"

using namespace Foam::constant::mathematical;

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(LarsenBorgnakkeVariableHardSphere, 0);
    addToRunTimeSelectionTable
    (
        BinaryCollisionModel,
        LarsenBorgnakkeVariableHardSphere,
        dictionary
    );
};

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::LarsenBorgnakkeVariableHardSphere::LarsenBorgnakkeVariableHardSphere
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
    vibrationalRelaxationCollisionNumber_
    (
        coeffDict_.lookupOrDefault<scalar>
        (
            "vibrationalRelaxationCollisionNumber", 
            0.0
        )
    ),
    electronicRelaxationCollisionNumber_
    (
        readScalar(coeffDict_.lookup("electronicRelaxationCollisionNumber"))
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::LarsenBorgnakkeVariableHardSphere::~LarsenBorgnakkeVariableHardSphere()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::LarsenBorgnakkeVariableHardSphere::active() const
{
    return true;
}


Foam::scalar Foam::LarsenBorgnakkeVariableHardSphere::sigmaTcR
(
    const dsmcParcel& pP,
    const dsmcParcel& pQ
) const
{
    const label typeIdP = pP.typeId();
    const label typeIdQ = pQ.typeId();

    const scalar dPQ =
        0.5
       *(
            cloud_.constProps(typeIdP).d()
          + cloud_.constProps(typeIdQ).d()
        );

    const scalar omegaPQ =
        0.5
       *(
            cloud_.constProps(typeIdP).omega()
          + cloud_.constProps(typeIdQ).omega()
        );

    const scalar cR = mag(pP.U() - pQ.U());

    if (cR < VSMALL)
    {
        return 0;
    }

    const scalar mP = cloud_.constProps(typeIdP).mass();
    const scalar mQ = cloud_.constProps(typeIdQ).mass();
    const scalar mR = mP*mQ/(mP + mQ);

    // calculating cross section = pi*dPQ^2, where dPQ is from Bird, eq. 4.79
    const scalar sigmaTPQ =
        pi*dPQ*dPQ
       *pow(2.0*physicoChemical::k.value()*Tref_/(mR*cR*cR), omegaPQ - 0.5)
       /exp(Foam::lgamma(2.5 - omegaPQ));

    return sigmaTPQ*cR;
}


void Foam::LarsenBorgnakkeVariableHardSphere::collide
(
    dsmcParcel& pP,
    dsmcParcel& pQ,
    label& cellI
)
{   
    Random& rndGen(cloud_.rndGen());
    
    const label typeIdP = pP.typeId();
    const label typeIdQ = pQ.typeId();
    
    vector& UP = pP.U();
    vector& UQ = pQ.U();
    
    const scalar mP = cloud_.constProps(typeIdP).mass();
    const scalar mQ = cloud_.constProps(typeIdQ).mass();
    const scalar mR = mP*mQ/(mP + mQ);
    const vector Ucm = (mP*UP + mQ*UQ)/(mP + mQ);
    const scalar cRsqr = magSqr(UP - UQ);

    //- Pre-collision relative translational energy
    scalar translationalEnergy = 0.5*mR*cRsqr;
    
    const scalar omegaPQ =
        0.5
       *(
            cloud_.constProps(typeIdP).omega()
          + cloud_.constProps(typeIdQ).omega()
        );
    
    relax(pP, translationalEnergy, omegaPQ);
    relax(pQ, translationalEnergy, omegaPQ);

    //- Rescale the translational energy
    const scalar cR = sqrt(2.0*translationalEnergy/mR);

    //- Variable Hard Sphere collision part (Eq. 3.14, p57)
    const scalar cosTheta = 2.0*rndGen.sample01<scalar>() - 1.0;
    const scalar sinTheta = sqrt(1.0 - cosTheta*cosTheta);
    const scalar phi = 2.0*pi*rndGen.sample01<scalar>();

    //- Post-collision relative velocity (Eq. 3.15, p57)
    const vector postCollisionRelU =
        cR
       *vector
        (
            cosTheta,
            sinTheta*cos(phi),
            sinTheta*sin(phi)
        );

    UP = Ucm + postCollisionRelU*mQ/(mP + mQ);

    UQ = Ucm - postCollisionRelU*mP/(mP + mQ);
    
    
    scalar collisionSeparation = sqrt(
            sqr(pP.position().x() - pQ.position().x()) +
            sqr(pP.position().y() - pQ.position().y())
    );
    
    cloud_.cellPropMeasurements().collisionSeparation()[cellI] += collisionSeparation;
    cloud_.cellPropMeasurements().nColls()[cellI]++;
    
    const label classificationP = pP.classification();
    const label classificationQ = pQ.classification();
    
    //- Class I molecule changes to class
    //  III molecule when it collides with either class II or class III
    //  molecules.
    
    if (classificationP == 0 && classificationQ == 1)
    {
        pP.classification() = 2;
    }
    
    if (classificationQ == 0 && classificationP == 1)
    {
        pQ.classification() = 2;
    }
    
    if (classificationP == 0 && classificationQ == 2)
    {
        pP.classification() = 2;
    }
    
    if (classificationQ == 0 && classificationP == 2)
    {
        pQ.classification() = 2;
    }
}


void Foam::LarsenBorgnakkeVariableHardSphere::relax
(
    dsmcParcel& p,
    scalar& translationalEnergy,
    const scalar omegaPQ,
    const bool postReaction
)
{
    const label typeIdP = p.typeId();
    const dsmcParcel::constantProperties& cP = cloud_.constProps(typeIdP);
    
    if (cP.type() == 0)
    {
        //- The particle is an electron, no energy to redistribute
        return void();
    }
    
    const scalar inverseRotationalCollisionNumber = 1.0/rotationalRelaxationCollisionNumber_;
    const scalar inverseElectronicCollisionNumber = 1.0/electronicRelaxationCollisionNumber_;
    
    scalar& ERotP = p.ERot();
    label& ELevelP = p.ELevel();
    
    //- Electronic energy mode for P
    if (inverseElectronicCollisionNumber > cloud_.rndGen().sample01<scalar>())
    { 
        const scalar preCollisionEEleP = cP.electronicEnergyList()[ELevelP];
        const label jMaxP = cP.nElectronicLevels();    
        const scalarList& EElistP = cP.electronicEnergyList();    
        const labelList& gListP = cP.electronicDegeneracyList();   
    
        //- Collision energy of particle P: relative translational energy 
        //   + pre-collision electronic energy
        const scalar EcP = translationalEnergy + preCollisionEEleP;
        
        const label postCollisionELevel = 
            cloud_.postCollisionElectronicEnergyLevel
                (
                    EcP,
                    jMaxP,
                    omegaPQ,
                    EElistP,
                    gListP
                );
                        
        ELevelP = postCollisionELevel;
        
        //- Relative translational energy after electronic energy exchange
        translationalEnergy = EcP - EElistP[ELevelP];
    }
            
    //- Vibrational energy mode for P
    if (cP.nVibrationalModes() > 0)
    {
        const scalarList& thetaVP = cP.thetaV();  
        const scalarList& thetaDP = cP.thetaD();
        const scalarList& ZrefP = cP.Zref();
        const scalarList& refTempZvP = cP.TrefZv();
        const scalarList& preCollisionEVibP = cP.eVib(p.vibLevel());
        
        forAll(thetaVP, i)
        {
            //- Collision energy of particle P: relative translational energy 
            //    + pre-collision vibrational energy
            const scalar EcP = translationalEnergy + preCollisionEVibP[i]; 

            //- Maximum possible quantum level (equation 3, Bird 2010)
            const label iMaxP = EcP/(physicoChemical::k.value()*thetaVP[i]); 

            if (iMaxP > 0)
            {       
                p.vibLevel()[i] = 
                    cloud_.postCollisionVibrationalEnergyLevel
                    (
                        postReaction,
                        p.vibLevel()[i],
                        iMaxP,
                        thetaVP[i],
                        thetaDP[i],
                        refTempZvP[i],
                        omegaPQ,
                        ZrefP[i],
                        EcP,
                        vibrationalRelaxationCollisionNumber_
                    );
                        
                translationalEnergy = EcP - cP.eVib_m(i, p.vibLevel()[i]);
            }
        }
    }
    
    //- Rotational energy mode for P
    const scalar rotationalDofP = cP.rotationalDegreesOfFreedom();
        
    // Larsen Borgnakke rotational energy redistribution part. Using the serial
    // application of the LB method, as per the INELRS subroutine in Bird's
    // DSMC0R.FOR
    if (rotationalDofP > VSMALL)
    {
         /*scalar particleProbabilityP = 
             ((zeta_T + 2.0*rotationalDofP)/(2.0*rotationalDofP))
             *(
                 1.0 - sqrt(
                             1.0 - (rotationalDofP/zeta_T)
                             *((zeta_T+rotationalDofP)/(zeta_T+2.0*rotationalDofP))
                             *(4.0/rotationalRelaxationCollisionNumber_)
                           )
              );
            
         Info << "particleProbabilityP = " << particleProbabilityP << endl;*/
       //if (particleProbabilityP > cloud_.rndGen().sample01<scalar>())  
        
        const scalar preCollisionERotP = ERotP;
        
        if (inverseRotationalCollisionNumber > cloud_.rndGen().sample01<scalar>())
        {
            const scalar EcP = translationalEnergy + preCollisionERotP;
            const scalar ChiB = 2.5 - omegaPQ;
            
            const scalar energyRatio = 
                cloud_.postCollisionRotationalEnergy(rotationalDofP, ChiB);

            ERotP = energyRatio*EcP;
        
            translationalEnergy = EcP - ERotP;
        }
    }
}


const Foam::dictionary&
Foam::LarsenBorgnakkeVariableHardSphere::coeffDict() const
{
    return coeffDict_;
}


// ************************************************************************* //
