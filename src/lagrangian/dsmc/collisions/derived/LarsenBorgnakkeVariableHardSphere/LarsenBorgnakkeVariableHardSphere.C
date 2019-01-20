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
    VariableHardSphere(dict, cloud),
    coeffDictLB_
    (
        dict.isDict(typeName + "Coeffs")
        ? dict.subDict(typeName + "Coeffs")
        : dictionary()
    ),
    rotationalRelaxationCollisionNumber_
    (
        coeffDictLB_.lookupOrDefault<scalar>
        (
            "rotationalRelaxationCollisionNumber",
            5.0
        )
    ),
    vibrationalRelaxationCollisionNumber_
    (
        coeffDictLB_.lookupOrDefault<scalar>
        (
            "vibrationalRelaxationCollisionNumber", 
            0.0
        )
    ),
    invZvFormulation_(2),
    electronicRelaxationCollisionNumber_
    (
        coeffDictLB_.lookupOrDefault<scalar>
        (
            "electronicRelaxationCollisionNumber",
            500.0
        )
    )
{
    const word inverseZvFormulationVersion =
        coeffDictLB_.lookupOrDefault<word>
        (
            "inverseZvFormulation", 
            word::null
        );
        
    if (inverseZvFormulationVersion == "pre-2008")
    {
        invZvFormulation_ = 0;
    }
    else if (inverseZvFormulationVersion == "2008")
    {
        invZvFormulation_ = 1;
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::LarsenBorgnakkeVariableHardSphere::~LarsenBorgnakkeVariableHardSphere()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::LarsenBorgnakkeVariableHardSphere::collide
(
    dsmcParcel& pP,
    dsmcParcel& pQ,
    const label cellI,
    scalar cR
)
{   
    const label typeIdP = pP.typeId();
    const label typeIdQ = pQ.typeId();
    
    vector& UP = pP.U();
    vector& UQ = pQ.U();
    
    const scalar mP = cloud_.constProps(typeIdP).mass();
    const scalar mQ = cloud_.constProps(typeIdQ).mass();
    const scalar mR = mP*mQ/(mP + mQ);
    
    const scalar cRsqr = magSqr(UP - UQ);

    //- Pre-collision relative translational energy
    scalar translationalEnergy = 0.5*mR*cRsqr;
    
    const scalar omegaPQ =
        0.5
       *(
            cloud_.constProps(typeIdP).omega()
          + cloud_.constProps(typeIdQ).omega()
        );
    
    redistribute(pP, translationalEnergy, omegaPQ);
    redistribute(pQ, translationalEnergy, omegaPQ);

    //- Rescale the translational energy
    cR = sqrt(2.0*translationalEnergy/mR);
    
    VariableHardSphere::scatter(pP, pQ, cellI, cR);
}


void Foam::LarsenBorgnakkeVariableHardSphere::redistribute
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
    
    const scalar inverseRotationalCollisionNumber =
        1.0/rotationalRelaxationCollisionNumber_;
    const scalar inverseElectronicCollisionNumber =
        1.0/electronicRelaxationCollisionNumber_;
    
    scalar& ERotP = p.ERot();
    label& ELevelP = p.ELevel();
    
    //- Electronic energy mode for P
    if (inverseElectronicCollisionNumber > cloud_.rndGen().sample01<scalar>())
    { 
        const label jMaxP = cP.nElectronicLevels();    
        const scalarList& EElistP = cP.electronicEnergyList();    
        const labelList& gListP = cP.electronicDegeneracyList(); 
        const scalar preCollisionEEleP = EElistP[ELevelP];
    
        //- Collision energy of particle P: relative translational energy 
        //   + pre-collision electronic energy
        const scalar EcP = translationalEnergy + preCollisionEEleP;
        
        ELevelP = 
            cloud_.postCollisionElectronicEnergyLevel
            (
                EcP,
                jMaxP,
                omegaPQ,
                EElistP,
                gListP
            );
                        
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
                        vibrationalRelaxationCollisionNumber_,
                        invZvFormulation_,
                        p.cell()
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
    if (rotationalDofP > 0)
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
    return coeffDictLB_;
}

// ************************************************************************* //
