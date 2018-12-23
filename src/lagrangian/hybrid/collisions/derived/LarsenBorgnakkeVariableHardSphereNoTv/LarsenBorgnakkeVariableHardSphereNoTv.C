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

#include "LarsenBorgnakkeVariableHardSphereNoTv.H"
#include "constants.H"
#include "addToRunTimeSelectionTable.H"

using namespace Foam::constant::mathematical;

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(LarsenBorgnakkeVariableHardSphereNoTv, 0);
    addToRunTimeSelectionTable(BinaryCollisionModel, LarsenBorgnakkeVariableHardSphereNoTv, dictionary);
};

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

//template <class CloudType>
Foam::LarsenBorgnakkeVariableHardSphereNoTv::LarsenBorgnakkeVariableHardSphereNoTv
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
Foam::LarsenBorgnakkeVariableHardSphereNoTv::~LarsenBorgnakkeVariableHardSphereNoTv()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::LarsenBorgnakkeVariableHardSphereNoTv::active() const
{
    return true;
}

//template <class CloudType>
Foam::scalar Foam::LarsenBorgnakkeVariableHardSphereNoTv::sigmaTcR
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
void Foam::LarsenBorgnakkeVariableHardSphereNoTv::collide
(
    dsmcParcel& pP,
    dsmcParcel& pQ
)
{   
    label typeIdP = pP.typeId();
    label typeIdQ = pQ.typeId();
    vector& UP = pP.U();
    vector& UQ = pQ.U();
    scalar& ERotP = pP.ERot();
    scalar& ERotQ = pQ.ERot();
    scalar& EVibP = pP.EVib();
    scalar& EVibQ = pQ.EVib();

    Random& rndGen(cloud_.rndGen());
    
    
  //   VIBRATIONAL ENERGY EXCHANGE - QUANTUM-KINETIC MODEL
    
    scalar preCollisionERotP = ERotP;

    scalar preCollisionERotQ = ERotQ;
    
    scalar preCollisionEVibP = EVibP;

    scalar preCollisionEVibQ = EVibQ;

    scalar rotationalDofP = cloud_.constProps(typeIdP).rotationalDegreesOfFreedom();

    scalar rotationalDofQ = cloud_.constProps(typeIdQ).rotationalDegreesOfFreedom();
    
    scalar vibrationalDofP = cloud_.constProps(typeIdP).nVibrationalModes();

    scalar vibrationalDofQ = cloud_.constProps(typeIdQ).nVibrationalModes();
    
    scalar thetaVP = cloud_.constProps(typeIdP).thetaV();
    
    scalar thetaVQ = cloud_.constProps(typeIdQ).thetaV();
    
    scalar thetaDP = cloud_.constProps(typeIdP).thetaD();
    
    scalar thetaDQ = cloud_.constProps(typeIdQ).thetaD();
    
    scalar ZrefP = cloud_.constProps(typeIdP).Zref();
    
    scalar ZrefQ = cloud_.constProps(typeIdQ).Zref();
    
    scalar refTempZvP = cloud_.constProps(typeIdP).TrefZv();
    
    scalar refTempZvQ = cloud_.constProps(typeIdQ).TrefZv();

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
    
//     scalar alphaPQ = 2.0/(omegaPQ-0.5);
    
//     scalar zeta_T = 4.0 - (4.0/alphaPQ);
    
//     Info << "alphaPQ = " << alphaPQ << endl;
//     Info << "zeta_T = " << zeta_T << endl;
            
    // Larsen Borgnakke rotational energy redistribution part.  Using the serial
    // application of the LB method, as per the INELRS subroutine in Bird's
    // DSMC0R.FOR
            
    if(vibrationalDofP > VSMALL)
    {
        // collision energy of particle P = relative translational energy + pre-collision vibrational energy
        scalar EcP = translationalEnergy + preCollisionEVibP; 

        // - maximum possible quantum level (equation 3, Bird 2010)
        label iMaxP = (EcP / (physicoChemical::k.value()*thetaVP)); 

        if(iMaxP > SMALL)
        {
            // - "quantised collision temperature" (equation 3, Bird 2010), denominator from Bird 5.42

            scalar TCollP = (iMaxP*thetaVP) / (3.5 - omegaPQ);
            
            scalar pow1 = pow((thetaDP/TCollP),0.33333) - 1.0;

            scalar pow2 = pow ((thetaDP/refTempZvP),0.33333) -1.0;

            // - vibrational collision number (equation 2, Bird 2010)
            scalar ZvP1 = pow((thetaDP/TCollP),omegaPQ); 
            
            scalar ZvP2 = pow(ZrefP*(pow((thetaDP/refTempZvP),(-1.0*omegaPQ))),(pow1/pow2));
            
//            scalar ZvP = ZvP1*ZvP2;
              scalar ZvP = 1.0;
            
            scalar inverseVibrationalCollisionNumberP = 1.0/ZvP;
        
            if(inverseVibrationalCollisionNumberP > rndGen.sample01<scalar>())
            {
                // post-collision quantum number
                label iDashP = 0; 
                scalar func = 0.0;
        
                do // acceptance - rejection 
                {
                    iDashP = rndGen.integer(0,iMaxP);
                    EVibP = iDashP*physicoChemical::k.value()*thetaVP;
                    
                    // - equation 5.61, Bird
                    func = pow((1.0 - (EVibP / EcP)),(1.5 - omegaPQ));

                } while( !(func > rndGen.sample01<scalar>()) );

                // relative translational energy after vibrational exchange
                translationalEnergy = EcP - EVibP;
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
        
        if (inverseRotationalCollisionNumber /*particleProbabilityP*/ > rndGen.sample01<scalar>())
        {
            scalar EcP = translationalEnergy + preCollisionERotP;
            
            scalar energyRatio = 0.0;
            
            if(rotationalDofP == 2.0)
            {
                energyRatio = 1.0 - pow(rndGen.sample01<scalar>(),(1.0/ChiB));
            }
            else
            {
                scalar ChiA = 0.5*rotationalDofP;
                
                energyRatio = cloud_.energyRatio(ChiA, ChiB);
            }

            ERotP = energyRatio*EcP;
        
            translationalEnergy = EcP - ERotP;
        }
    }   

      
    if(vibrationalDofQ > VSMALL)
    {
        scalar EcQ = translationalEnergy + preCollisionEVibQ;
        label iMaxQ = (EcQ /(physicoChemical::k.value()*thetaVQ));

        if(iMaxQ > SMALL)
        {
            // - Bird equations 5.42 and 11.34 gave this denominator
            scalar TCollQ = (iMaxQ*thetaVQ) / (3.5 - omegaPQ); 
            
            scalar pow1 = pow((thetaDQ/TCollQ),0.333) - 1.0;

            scalar pow2 = pow ((thetaDQ/refTempZvQ),0.333) -1.0;
            
            // - vibrational collision number (equation 2, Bird 2010)
            scalar ZvQ1 = pow((thetaDQ/TCollQ),omegaPQ); 
            
            scalar ZvQ2 = pow(ZrefQ*(pow((thetaDQ/refTempZvQ),(-1.0*omegaPQ))),(pow1/pow2));
            
//            scalar ZvQ = ZvQ1*ZvQ2;
             scalar ZvQ = 1.0;
                
            scalar inverseVibrationalCollisionNumberQ = 1.0/ZvQ;
        
            if(inverseVibrationalCollisionNumberQ > rndGen.sample01<scalar>())
            {
                label iDashQ = 0; // post-collision quantum number
                scalar func = 0.0;

                do // acceptance - rejection 
                {
                    iDashQ = rndGen.integer(0,iMaxQ);
                    EVibQ = iDashQ*physicoChemical::k.value()*thetaVQ;
                    func = pow((1.0 - (EVibQ / EcQ)),(1.5 - omegaPQ));
            
                } while( !(func > rndGen.sample01<scalar>()) );
        
                translationalEnergy = EcQ - EVibQ;
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
        
        if (inverseRotationalCollisionNumber /*particleProbabilityQ*/ > rndGen.sample01<scalar>())
        {
            scalar EcQ = translationalEnergy + preCollisionERotQ;
            
            scalar energyRatio = 0.0;
            
            if(rotationalDofQ == 2.0)
            {
                energyRatio = 1.0 - pow(rndGen.sample01<scalar>(),(1.0/ChiB));   
            }
            else
            {
                scalar ChiA = 0.5*rotationalDofQ;
                
                energyRatio = cloud_.energyRatio(ChiA, ChiB);
            }

            ERotQ = energyRatio*EcQ;
        
            translationalEnergy = EcQ - ERotQ;
        }
    }

    // Rescale the translational energy
    scalar cR = sqrt((2.0*translationalEnergy)/mR);

    // Variable Hard Sphere collision part

    scalar cosTheta = 2.0*rndGen.sample01<scalar>() - 1.0;

    scalar sinTheta = sqrt(1.0 - cosTheta*cosTheta);

    scalar phi = 2.0*pi*rndGen.sample01<scalar>();

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
}

const Foam::dictionary& Foam::LarsenBorgnakkeVariableHardSphereNoTv::coeffDict() const
{
    return coeffDict_;
}


// ************************************************************************* //
