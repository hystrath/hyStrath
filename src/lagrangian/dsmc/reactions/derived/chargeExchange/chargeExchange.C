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

#include "chargeExchange.H"
#include "addToRunTimeSelectionTable.H"
#include "fvc.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(chargeExchange, 0);

addToRunTimeSelectionTable(dsmcReaction, chargeExchange, dictionary);



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
chargeExchange::chargeExchange
(
    Time& t,
    dsmcCloud& cloud,
    const dictionary& dict
)
:
    dsmcReaction(t, cloud, dict),
    propsDict_(dict.subDict(typeName + "Properties")),
    reactantIds_(),
    chargeExchangeProductIds_(),
    dissociationProductIds_(),
    ionisationProductIds_(),
    reactionName_(propsDict_.lookup("reactionName")),
    dissociationPossible_(Switch(propsDict_.lookup("dissociationPossible"))),
    heatOfReactionDissoc_(0.0),
    heatOfReactionIon_(readScalar(propsDict_.lookup("heatOfReactionIonisation"))),
    heatOfReactionExch_(readScalar(propsDict_.lookup("heatOfReactionExchange"))),
    activationEnergy_(0.0),
    aCoeff_(readScalar(propsDict_.lookup("aCoeff"))),
    bCoeff_(readScalar(propsDict_.lookup("bCoeff"))),
    relax_(true),
    allowSplitting_(false),
    writeRatesToTerminal_(false),
    nChargeExchangeReactions_(0),
    nIonisationReactions_(0),
    nDissociationReactions_(0),
    nChargeExchangeReactionsPerTimeStep_(0),
    nDissociationReactionsPerTimeStep_(0),
    nIonisationReactionsPerTimeStep_(0),
    volume_(0.0),
    numberDensities_(2, 0.0)
{

}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

chargeExchange::~chargeExchange()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void chargeExchange::initialConfiguration()
{
    setProperties();
}

void chargeExchange::setProperties()
{
    // reading in reactants

    const List<word> reactantMolecules (propsDict_.lookup("reactants"));

    if(reactantMolecules.size() != 2)
    {
        FatalErrorIn("chargeExchange::setProperties()")
            << "There should be two reactants atoms, instead of " 
            << reactantMolecules.size() << nl 
            << exit(FatalError);
    }

    reactantIds_.setSize(reactantMolecules.size(), -1);

    allowSplitting_ = Switch(propsDict_.lookup("allowSplitting"));
    
    writeRatesToTerminal_ = Switch(propsDict_.lookup("writeRatesToTerminal"));

    forAll(reactantIds_, r)
    {
        reactantIds_[r] = findIndex(cloud_.typeIdList(), reactantMolecules[r]);

        // check that reactants belong to the typeIdList (constant/dsmcProperties)
        if(reactantIds_[r] == -1)
        {
            FatalErrorIn("chargeExchange::setProperties()")
                << "Cannot find type id: " << reactantMolecules[r] << nl 
                << exit(FatalError);
        }

        // check that reactants are 'ATOMS' or 'MOLECULES', not 'ELECTRONS'

        const label& charge = cloud_.constProps(reactantIds_[r]).charge();
    
        if(charge == -1)
        {
            FatalErrorIn("chargeExchange::setProperties()")
                << "Reactant must not be an electron: " << reactantMolecules[r] 
                << nl 
                << exit(FatalError);
        }
        
        const label& vDof = cloud_.constProps(reactantIds_[r]).vibrationalDegreesOfFreedom();

        if(vDof > 1)
        {
             FatalErrorIn("chargeExchange::setProperties()")
            << "Reactions are currently only implemented for monatomic and diatomic species"
            << " This is a polyatomic:" << reactantMolecules[r] 
            << nl 
            << exit(FatalError);
        }
    }

    // reading in charge exchange products

    const List<word> chargeExchangeProductMolecules (propsDict_.lookup("chargeExchangeProducts"));
    
    chargeExchangeProductIds_.setSize(chargeExchangeProductMolecules.size(), -1);
    
    forAll(chargeExchangeProductIds_, i)
    {
        chargeExchangeProductIds_[i] = findIndex(cloud_.typeIdList(), chargeExchangeProductMolecules[i]);
        
        // check that products belong to the typeIdList (constant/dsmcProperties)
        if(chargeExchangeProductIds_[i] == -1)
        {
            FatalErrorIn("chargeExchange::setProperties()")
                << "Cannot find type id: " << chargeExchangeProductMolecules[i] << nl 
                << exit(FatalError);
        }
        
        // check that products are a 'MOLECULE' or an 'ATOM'
        
        const label& charge = cloud_.constProps(chargeExchangeProductIds_[i]).charge();

        if(charge == -1)
        {
            FatalErrorIn("chargeExchange::setProperties()")
                << "Products cannot be an electron: " << chargeExchangeProductMolecules 
                << nl 
                << exit(FatalError);
        }
    }
    
    // reading in dissociation products

    if(dissociationPossible_)
    {
        const List<word> dissociationProductMolecules (propsDict_.lookup("dissociationProducts"));
        
        dissociationProductIds_.setSize(dissociationProductMolecules.size(), -1);
        
        forAll(dissociationProductIds_, i)
        {
            dissociationProductIds_[i] = findIndex(cloud_.typeIdList(), dissociationProductMolecules[i]);
            
            // check that products belong to the typeIdList (constant/dsmcProperties)
            if(dissociationProductIds_[i] == -1)
            {
                FatalErrorIn("chargeExchange::setProperties()")
                    << "Cannot find type id: " << dissociationProductMolecules[i] << nl 
                    << exit(FatalError);
            }
            
            // check that products are ATOMS
            
            const scalar& rDof = cloud_.constProps(dissociationProductIds_[i]).rotationalDegreesOfFreedom();

            if(rDof > VSMALL)
            {
                FatalErrorIn("chargeExchange::setProperties()")
//                     << "Dissociation products must be atoms: " << dissociationProductMolecules 
                    << nl 
                    << exit(FatalError);
            }
        }
    }
    
    // reading in ionisation products

    const List<word> ionisationProductMolecules (propsDict_.lookup("ionisationProducts"));
    
    ionisationProductIds_.setSize(ionisationProductMolecules.size(), -1);
    
    forAll(ionisationProductIds_, i)
    {
        ionisationProductIds_[i] = findIndex(cloud_.typeIdList(), ionisationProductMolecules[i]);
        
        // check that products belong to the typeIdList (constant/dsmcProperties)
        if(ionisationProductIds_[i] == -1)
        {
            FatalErrorIn("chargeExchange::setProperties()")
                << "Cannot find type id: " << ionisationProductMolecules[i] << nl 
                << exit(FatalError);
        }
    }
    
    // check that second product is an electron
        
    const label& charge = cloud_.constProps(ionisationProductIds_[1]).charge();

    if(charge != -1)
    {
        FatalErrorIn("chargeExchange::setProperties()")
            << "Second ionisation product must be an electron: " << ionisationProductMolecules 
            << nl 
            << exit(FatalError);
    }
    
    activationEnergy_ *= physicoChemical::k.value();
    
    if(dissociationPossible_)
    {
        heatOfReactionDissoc_ = readScalar(propsDict_.lookup("heatOfReactionDissociation"));
    }

}

bool chargeExchange::tryReactMolecules(const label& typeIdP, const label& typeIdQ)
{
    label reactantPId = findIndex(reactantIds_, typeIdP);
    label reactantQId = findIndex(reactantIds_, typeIdQ);

    if(reactantPId != reactantQId)
    {
        if
        (
            (reactantPId != -1) &&
            (reactantQId != -1) 
        )
        {
            return true;
        }
    }

    if
    (
        (reactantPId != -1) &&
        (reactantQId != -1) &&
        (reactantPId != reactantQId)
    )
    {
        return true;
    }
    else
    {
        return false;
    }
}

void chargeExchange::reaction
(
    dsmcParcel& p,
    dsmcParcel& q,
    const DynamicList<label>& candidateList,
    const List<DynamicList<label> >& candidateSubList,
    const label& candidateP,
    const List<label>& whichSubCell
)
{
}


void chargeExchange::reaction
(
    dsmcParcel& p,
    dsmcParcel& q
)
{
    label typeIdP = p.typeId();
    label typeIdQ = q.typeId();
    
    if(typeIdP == reactantIds_[0] &&  typeIdQ == reactantIds_[1])
    {
        relax_ = true;
        
        scalar totalReactionProbability = 0.0;
        scalarList reactionProbabilities(3, 0.0);
        
        vector UP = p.U();
        vector UQ = q.U();
        
        scalarList EElistP = cloud_.constProps(typeIdP).electronicEnergyList();
        labelList gListP = cloud_.constProps(typeIdP).degeneracyList();
        
        scalarList EElistQ = cloud_.constProps(typeIdQ).electronicEnergyList();
        labelList gListQ = cloud_.constProps(typeIdQ).degeneracyList();
        
        scalar ERotP = p.ERot();
        scalar EVibP = 0.0;
        
        if(cloud_.constProps(typeIdP).vibrationalDegreesOfFreedom() > 0)
        {
            EVibP = 
                p.vibLevel()[0]
                *cloud_.constProps(typeIdP).thetaV()[0]
                *physicoChemical::k.value();
        }
        scalar EElecP = EElistP[p.ELevel()];
        
        scalar ERotQ = q.ERot();
        scalar EVibQ = 0.0;
        
        if(cloud_.constProps(typeIdQ).vibrationalDegreesOfFreedom() > 0)
        {
            EVibQ = 
                q.vibLevel()[0]
                *cloud_.constProps(typeIdQ).thetaV()[0]
                *physicoChemical::k.value();
        }
        
        scalar EElecQ = EElistQ[q.ELevel()];
        
        label maxElectronicLevelP = cloud_.constProps(typeIdP).numberOfElectronicLevels();

        scalar mP = cloud_.constProps(typeIdP).mass();
        scalar mQ = cloud_.constProps(typeIdQ).mass();

        scalar mR = mP*mQ/(mP + mQ);
        scalar cRsqr = magSqr(UP - UQ);
        scalar translationalEnergy = 0.5*mR*cRsqr;
        
        scalar omegaPQ =
        0.5
        *(
            cloud_.constProps(typeIdP).omega()
            + cloud_.constProps(typeIdQ).omega()
        );
        
        scalar ChiB = 2.5 - omegaPQ;
        
        bool dissocReactionP = false;
        bool ionisationReactionP = false;
        bool chargeExchange = false;
                
        //3 reactions possible
        // 1. Dissociation of P
        // 2. Ionisation of P
        // 3. Charge exchange

        
        scalar EcP = 0.0;
        scalar TColl = 0.0;
        //check for dissociation
        
        if(dissociationPossible_)
        {
            if(cloud_.constProps(typeIdP).rotationalDegreesOfFreedom() > 0)
            {
                label idP = 
                    cloud_.constProps(typeIdP).charDissQuantumLevel()[0];
                label imaxP = 0;
                
                //Only diatomics can be dissociated
                EcP = translationalEnergy + EVibP;

                imaxP = EcP/(physicoChemical::k.value()*cloud_.constProps(typeIdP).thetaV()[0]);
                
                if(imaxP-idP > 0)
                {
                    //Dissociation can occur
                    totalReactionProbability += 1.0;
                    reactionProbabilities[0] = 1.0;
                }
            }
        }
        
        scalar ionisationEnergy = cloud_.constProps(typeIdP).ionisationTemperature()*physicoChemical::k.value();
        
        // calculate if an ionisation of species P is possible
        EcP = translationalEnergy + EElecP;

        if((EcP - ionisationEnergy) > VSMALL)
        {
            //Ionisation can occur
            totalReactionProbability += 1.0;
            reactionProbabilities[1] = 1.0;
        }
        
        scalar heatOfReactionExchJoules = heatOfReactionExch_*physicoChemical::k.value();
        
        
        
        EcP = translationalEnergy + EElecP;
        
        scalar aDash = aCoeff_*(pow(2.5 - omegaPQ, bCoeff_)*exp(lgamma(2.5 - omegaPQ))/exp(lgamma(2.5 - omegaPQ + bCoeff_)));
        
        TColl = (translationalEnergy/(physicoChemical::k.value()))/(2.5 - omegaPQ);
        
        scalar activationEnergy = activationEnergy_ + (aDash*pow((TColl/273.0) , bCoeff_) * fabs(heatOfReactionExchJoules));
        
//         scalar activationEnergy = activationEnergy_ + (aCoeff_*pow((5000/273.0) , bCoeff_) * fabs(heatOfReactionExchJoules));
        
        if(heatOfReactionExchJoules < 0.0)
        {
            activationEnergy -= heatOfReactionExchJoules;
        }

        if(EcP > activationEnergy)
        {
            label keyElectronicLevel = -1;
            
            for(label i = 0 ; i < maxElectronicLevelP ; i++)
            {           
                scalar electronicEnergy = EElistP[i];
                
                if(electronicEnergy > activationEnergy)
                {
                    break;
                }

                keyElectronicLevel++;
            }
            
            EcP = translationalEnergy + EElecP /*+ heatOfReactionExchJoules*/;
            
            label trialELevel = cloud_.postCollisionElectronicEnergyLevel
                            (
                                EcP,
                                maxElectronicLevelP,
                                omegaPQ,
                                EElistP,
                                gListP
                            );
                            
            if(trialELevel == keyElectronicLevel)
            {
                scalar prob = 0.0;
                
                label nPossStates = 0;
                    
                if(maxElectronicLevelP == 1)
                {
                    nPossStates = gListP[0];
                }
                else
                {
                    forAll(EElistP, n)
                    {
                        if(EcP > EElistP[n])
                        {
                            nPossStates += gListP[n];
                        }
                    }
                }
                
                label nState = ceil(cloud_.rndGen().scalar01()*(nPossStates));
                label nAvailableStates = 0;
                label nLevel = -1;
                
                forAll(EElistP, n)
                {
                    nAvailableStates += gListP[n];
                    
                    if(nState <= nAvailableStates && nLevel < 0)
                    {
                        nLevel = n;
                    }
                }
                
                //Calculate the probability of it occuring
                scalar summation = 0.0;

                for(label i = 0 ; i <= nLevel; i++)
                { 
                    summation += gListP[i]*pow((EcP - EElistP[i]),1.5-omegaPQ);
                }
                
                prob = (gListP[trialELevel]*pow((EcP - EElistP[trialELevel]),1.5-omegaPQ))/summation;
                
                if(prob > cloud_.rndGen().scalar01())
                {
                    //Charge exchange can occur
                    totalReactionProbability += prob;
                    reactionProbabilities[2] = prob;
                }
            }
        }
        
        //Decide if a reaction is to occur
        
        if(totalReactionProbability > cloud_.rndGen().scalar01())
        {
            //A chemical reaction is to occur, choose which one
            
            scalarList normalisedProbabilities(reactionProbabilities.size(), 0.0);
            scalar cumulativeProbability = 0.0;
            
            normalisedProbabilities = reactionProbabilities/totalReactionProbability;
            
            forAll(normalisedProbabilities, i)
            {              
                //If current reaction can't occur, don't check for it
                if(normalisedProbabilities[i] > VSMALL)
                {
                    cumulativeProbability += normalisedProbabilities[i];
                    
                    if(cumulativeProbability > cloud_.rndGen().scalar01())
                    {
                        //Current reaction is to occur
                        
                        if(i == 0)
                        {
                            //Dissociation is to occur
                            dissocReactionP = true;
                            break;
                        }
                        if(i == 1)
                        {
                            //Ionisation reaction is to occur
                            ionisationReactionP = true;
                            break;
                        }
                        if(i == 2)
                        {
                            //Dissociation is to occur
                            chargeExchange = true;
                            break;
                        }
                    }
                }
            }
        }
        
        if(dissocReactionP)
        {
            nDissociationReactions_++;
            nDissociationReactionsPerTimeStep_++;
            
            if(allowSplitting_)
            {
                relax_ = false;
                
                const scalar& heatOfReactionDissociationJoules = heatOfReactionDissoc_*physicoChemical::k.value();
                
                translationalEnergy = translationalEnergy + heatOfReactionDissociationJoules + EVibP;
                
                scalar thetaVQ = 0;
                scalar thetaDQ = 0;
                scalar ZrefQ = 0;
                scalar refTempZvQ = 0;
                
                if(cloud_.constProps(typeIdQ).vibrationalDegreesOfFreedom() > 0)
                {
                    thetaVQ = cloud_.constProps(typeIdQ).thetaV()[0];
                    thetaDQ = cloud_.constProps(typeIdQ).thetaD()[0];
                    ZrefQ = cloud_.constProps(typeIdQ).Zref()[0];
                    refTempZvQ = cloud_.constProps(typeIdQ).TrefZv()[0];
                }
                
                scalar jMaxQ = cloud_.constProps(typeIdQ).numberOfElectronicLevels()-1;
                scalar rotationalDofQ = cloud_.constProps(typeIdQ).rotationalDegreesOfFreedom();
                
                
                translationalEnergy += EElecQ;
                    
                label ELevelQ = cloud_.postCollisionElectronicEnergyLevel
                                (
                                    translationalEnergy,
                                    jMaxQ,
                                    omegaPQ,
                                    EElistQ,
                                    gListQ
                                );
                                
                translationalEnergy -= EElistQ[ELevelQ];
                
                label vibLevelQ = 0;

                if(rotationalDofQ > VSMALL)
                {
                    translationalEnergy += EVibQ;
                    
                    label iMax = (translationalEnergy / (physicoChemical::k.value()*thetaVQ));
                    
                    vibLevelQ = cloud_.postCollisionVibrationalEnergyLevel
                                    (
                                            true,
                                            q.vibLevel()[0],
                                            iMax,
                                            thetaVQ,
                                            thetaDQ,
                                            refTempZvQ,
                                            omegaPQ,
                                            ZrefQ,
                                            translationalEnergy
                                        );
                                    
                    translationalEnergy -= vibLevelQ*thetaVQ*physicoChemical::k.value();
                                    
                    translationalEnergy += ERotQ;
                    
                    ERotQ = translationalEnergy*cloud_.postCollisionRotationalEnergy(rotationalDofQ,ChiB);
                            
                    translationalEnergy -= ERotQ;
                }
                
                scalar relVelNonDissoMol = sqrt(2.0*translationalEnergy/mR);

                //center of mass velocity of all particles
                vector Ucm = (mP*UP + mQ*UQ)/(mP + mQ);    

                // Variable Hard Sphere collision part

                scalar cosTheta = 2.0*cloud_.rndGen().scalar01() - 1.0;
            
                scalar sinTheta = sqrt(1.0 - cosTheta*cosTheta);
            
                scalar phi = twoPi*cloud_.rndGen().scalar01();
            
                vector postCollisionRelU =
                    relVelNonDissoMol
                    *vector
                        (
                            cosTheta,
                            sinTheta*cos(phi),
                            sinTheta*sin(phi)
                        );

                UP = Ucm + (postCollisionRelU*mQ/(mP + mQ)); 
                UQ = Ucm - (postCollisionRelU*mP/(mP + mQ));

                const label& typeId1 = dissociationProductIds_[0];
                const label& typeId2 = dissociationProductIds_[1];
                
                //Mass of Product one and two
                scalar mP1 = cloud_.constProps(typeId1).mass();
                scalar mP2 = cloud_.constProps(typeId2).mass();
                
                scalar mRatoms = mP1*mP2/(mP1 + mP2);
                
                translationalEnergy = ERotP + EElecP;
                
                scalar cRatoms = sqrt(2.0*translationalEnergy/mRatoms);

                // Variable Hard Sphere collision part
                scalar cosTheta2 = 2.0*cloud_.rndGen().scalar01() - 1.0;
            
                scalar sinTheta2 = sqrt(1.0 - cosTheta2*cosTheta2);
            
                scalar phi2 = twoPi*cloud_.rndGen().scalar01();
            
                vector postCollisionRelU2 = cRatoms
                *vector
                    (
                        cosTheta2,
                        sinTheta2*cos(phi2),
                        sinTheta2*sin(phi2)
                    );


                vector uP1 = UP + postCollisionRelU2*mP2/(mP1 + mP2);
                vector uP2 = UP - postCollisionRelU2*mP1/(mP1 + mP2);

                // Q remains NON-DISSOCIATED.
                q.U() = UQ;
                q.ERot() = ERotQ;
                if(rotationalDofQ > VSMALL)
                {
                    q.vibLevel().setSize(1, vibLevelQ);
                }
                else 
                {
                    q.vibLevel().setSize(0,0);
                }
                q.ELevel() = ELevelQ;

                // Molecule P will dissociate
                vector position = p.position();
                
                label cell = -1;
                label tetFace = -1;
                label tetPt = -1;

                mesh_.findCellFacePt
                (
                    position,
                    cell,
                    tetFace,
                    tetPt
                );
                
                p.typeId() = typeId1;
                p.U() = uP1;
                p.vibLevel().setSize(0,0);
                p.ERot() = 0.0;
                p.ELevel() = 0;
                
                label classificationP = p.classification();
                scalar RWF = p.RWF();
                labelList vibLevel;
                
                // insert new product 2
                cloud_.addNewParcel
                (
                    position,
                    uP2,
                    RWF,
                    0.0,
                    0,
                    cell,
                    tetFace,
                    tetPt,
                    typeId2,
                    -1,
                    classificationP,
                    vibLevel
                );
            }
        }
        
        if(ionisationReactionP)
        {
            nIonisationReactions_++;
            nIonisationReactionsPerTimeStep_++;
            
            if(allowSplitting_)
            {
                relax_ = false;
                
                const scalar& heatOfReactionIonisationJoules = heatOfReactionIon_*physicoChemical::k.value();
                
                translationalEnergy = translationalEnergy + heatOfReactionIonisationJoules + EElecP;
                
                scalar thetaVQ = 0;
                scalar thetaDQ = 0;
                scalar ZrefQ = 0;
                scalar refTempZvQ = 0;
                
                if(cloud_.constProps(typeIdQ).vibrationalDegreesOfFreedom() > 0)
                {
                    thetaVQ = cloud_.constProps(typeIdQ).thetaV()[0];
                    thetaDQ = cloud_.constProps(typeIdQ).thetaD()[0];
                    ZrefQ = cloud_.constProps(typeIdQ).Zref()[0];
                    refTempZvQ = cloud_.constProps(typeIdQ).TrefZv()[0];
                }
               
                scalar jMaxQ = cloud_.constProps(typeIdQ).numberOfElectronicLevels()-1;
                scalar rotationalDofP = cloud_.constProps(typeIdP).rotationalDegreesOfFreedom();
                scalar rotationalDofQ = cloud_.constProps(typeIdQ).rotationalDegreesOfFreedom();
 
                translationalEnergy += EElecQ;
                    
                label ELevelQ = cloud_.postCollisionElectronicEnergyLevel
                                (
                                    translationalEnergy,
                                    jMaxQ,
                                    omegaPQ,
                                    EElistQ,
                                    gListQ
                                );
                                
                translationalEnergy -= EElistQ[ELevelQ];
                
                label vibLevelQ = 0;

                if(rotationalDofQ > VSMALL)
                {
                    translationalEnergy += EVibQ;
                    
                    label iMax = (translationalEnergy / (physicoChemical::k.value()*thetaVQ));
                    
                    vibLevelQ = cloud_.postCollisionVibrationalEnergyLevel
                                    (
                                            true,
                                            q.vibLevel()[0],
                                            iMax,
                                            thetaVQ,
                                            thetaDQ,
                                            refTempZvQ,
                                            omegaPQ,
                                            ZrefQ,
                                            translationalEnergy
                                        );
                                    
                    translationalEnergy -= vibLevelQ*thetaVQ*physicoChemical::k.value();
                                    
                    translationalEnergy += ERotQ;
                    
                    ERotQ = translationalEnergy*cloud_.postCollisionRotationalEnergy(rotationalDofQ,ChiB);
                            
                    translationalEnergy -= ERotQ;
                }
                
                scalar relVelNonDissoMol = sqrt(2.0*translationalEnergy/mR);

                //center of mass velocity of all particles
                vector Ucm = (mP*UP + mQ*UQ)/(mP + mQ);    

                // Variable Hard Sphere collision part

                scalar cosTheta = 2.0*cloud_.rndGen().scalar01() - 1.0;
            
                scalar sinTheta = sqrt(1.0 - cosTheta*cosTheta);
            
                scalar phi = twoPi*cloud_.rndGen().scalar01();
            
                vector postCollisionRelU =
                    relVelNonDissoMol
                    *vector
                        (
                            cosTheta,
                            sinTheta*cos(phi),
                            sinTheta*sin(phi)
                        );

                UP = Ucm + (postCollisionRelU*mQ/(mP + mQ)); 
                UQ = Ucm - (postCollisionRelU*mP/(mP + mQ));

                const label& typeId1 = ionisationProductIds_[0];
                const label& typeId2 = ionisationProductIds_[1];
                
                //Mass of Product one and two
                scalar mP1 = cloud_.constProps(typeId1).mass();
                scalar mP2 = cloud_.constProps(typeId2).mass();
                
                scalar mRatoms = mP1*mP2/(mP1 + mP2);
                
                translationalEnergy = ERotP + EVibP;
                
                scalar cRatoms = sqrt(2.0*translationalEnergy/mRatoms);

                // Variable Hard Sphere collision part
                scalar cosTheta2 = 2.0*cloud_.rndGen().scalar01() - 1.0;
            
                scalar sinTheta2 = sqrt(1.0 - cosTheta2*cosTheta2);
            
                scalar phi2 = twoPi*cloud_.rndGen().scalar01();
            
                vector postCollisionRelU2 = cRatoms
                *vector
                    (
                        cosTheta2,
                        sinTheta2*cos(phi2),
                        sinTheta2*sin(phi2)
                    );


                vector uP1 = UP + postCollisionRelU2*mP2/(mP1 + mP2);
                vector uP2 = UP - postCollisionRelU2*mP1/(mP1 + mP2);

                // Q remains NON-DISSOCIATED.
                q.U() = UQ;
                q.ERot() = ERotQ;
                if(rotationalDofQ > VSMALL)
                {
                    q.vibLevel()[0] = vibLevelQ;
                }
                else
                {
                    q.vibLevel().setSize(0,0);
                }
                q.ELevel() = ELevelQ;

                // Molecule P will ionise.
                vector position = p.position();
                
                label cell = -1;
                label tetFace = -1;
                label tetPt = -1;

                mesh_.findCellFacePt
                (
                    position,
                    cell,
                    tetFace,
                    tetPt
                );
                
                p.typeId() = typeId1;
                p.U() = uP1;
                if(rotationalDofP > VSMALL)
                {
                    p.vibLevel().setSize(1,0);
                }       
                else
                {
                    p.vibLevel().setSize(0,0);
                } 
                p.ERot() = 0.0;
                p.ELevel() = 0;
                
                label classificationP = p.classification();
                scalar RWF = p.RWF();
                labelList vibLevel;
                
                // insert new product 2
                cloud_.addNewParcel
                (
                    position,
                    uP2,
                    RWF,
                    0.0,
                    0,
                    cell,
                    tetFace,
                    tetPt,
                    typeId2,
                    -1,
                    classificationP,
                    vibLevel
                );
            }
        }
        
        if(chargeExchange)
        {
            nChargeExchangeReactions_++;
            nChargeExchangeReactionsPerTimeStep_++;
            
            if(allowSplitting_)
            {
                relax_ = false;
                
                translationalEnergy = translationalEnergy + heatOfReactionExchJoules;
                
                translationalEnergy += ERotP + EVibP + EElecP + ERotQ + EVibQ + EElecQ;
                
                scalar relVel = sqrt((2.0*translationalEnergy)/mR);
                    
                // centre of mass velocity of molecules (pre-split)
                vector Ucm = (mP*UP + mQ*UQ)/(mP + mQ);

                // Variable Hard Sphere collision part for collision of molecules
                scalar cosTheta = 2.0*cloud_.rndGen().scalar01() - 1.0;
            
                scalar sinTheta = sqrt(1.0 - cosTheta*cosTheta);
            
                scalar phi = twoPi*cloud_.rndGen().scalar01();
            
                vector postCollisionRelU =
                    relVel
                    *vector
                        (
                            cosTheta,
                            sinTheta*cos(phi),
                            sinTheta*sin(phi)
                        );
                        
                mP = cloud_.constProps(chargeExchangeProductIds_[0]).mass();
                mQ = cloud_.constProps(chargeExchangeProductIds_[1]).mass();
        
                UP = Ucm + (postCollisionRelU*mQ/(mP + mQ));
                UQ = Ucm - (postCollisionRelU*mP/(mP + mQ));
                
                p.typeId() = chargeExchangeProductIds_[0];
                p.U() = UP;
                p.ERot() = 0.0;
                if(cloud_.constProps(chargeExchangeProductIds_[0]).rotationalDegreesOfFreedom() > VSMALL)
                {
                    p.vibLevel().setSize(1,0);
                }       
                else
                {
                    p.vibLevel().setSize(0,0);
                } 
                p.ELevel() = 0;
                
                q.typeId() = chargeExchangeProductIds_[1];
                q.U() = UQ;
                q.ERot() = 0.0;
                if(cloud_.constProps(chargeExchangeProductIds_[1]).rotationalDegreesOfFreedom() > VSMALL)
                {
                    q.vibLevel().setSize(1,0);
                }       
                else
                {
                    q.vibLevel().setSize(0,0);
                } 
                q.ELevel() = 0;
            }              
        }
    }
    
    if(typeIdP == reactantIds_[1] &&  typeIdQ == reactantIds_[0])
    {
        relax_ = true;
        
        scalar totalReactionProbability = 0.0;
        scalarList reactionProbabilities(3, 0.0);
        
        vector UP = p.U();
        vector UQ = q.U();
        
        scalarList EElistP = cloud_.constProps(typeIdP).electronicEnergyList();
        labelList gListP = cloud_.constProps(typeIdP).degeneracyList();
        
        scalarList EElistQ = cloud_.constProps(typeIdQ).electronicEnergyList();
        labelList gListQ = cloud_.constProps(typeIdQ).degeneracyList();
        
        scalar ERotP = p.ERot();
        scalar EVibP = 0.0;
        
        if(cloud_.constProps(typeIdP).vibrationalDegreesOfFreedom() > 0)
        {
            EVibP = 
                p.vibLevel()[0]
                *cloud_.constProps(typeIdP).thetaV()[0]
                *physicoChemical::k.value();
        }
        
        scalar EElecP = EElistP[p.ELevel()];
        
        scalar ERotQ = q.ERot();
        scalar EVibQ = 0.0;
        
        if(cloud_.constProps(typeIdQ).vibrationalDegreesOfFreedom() > 0)
        {
            EVibQ = 
                q.vibLevel()[0]
                *cloud_.constProps(typeIdQ).thetaV()[0]
                *physicoChemical::k.value();
        }
        scalar EElecQ = EElistQ[q.ELevel()];
        
        label maxElectronicLevelQ = cloud_.constProps(typeIdQ).numberOfElectronicLevels();

        scalar mP = cloud_.constProps(typeIdP).mass();
        scalar mQ = cloud_.constProps(typeIdQ).mass();

        scalar mR = mP*mQ/(mP + mQ);
        scalar cRsqr = magSqr(UP - UQ);
        scalar translationalEnergy = 0.5*mR*cRsqr;
        
        scalar omegaPQ =
        0.5
        *(
            cloud_.constProps(typeIdP).omega()
            + cloud_.constProps(typeIdQ).omega()
        );
        
        scalar ChiB = 2.5 - omegaPQ;
        
        bool dissocReactionQ = false;
        bool ionisationReactionQ = false;
        bool chargeExchange = false;
                
        //3 reactions possible
        // 1. Dissociation of P
        // 2. Ionisation of P
        // 3. Charge exchange

        
        scalar EcQ = 0.0;
        scalar TColl = 0.0;
        
        //check for dissociation
        
        if(dissociationPossible_)
        {
            if(cloud_.constProps(typeIdQ).rotationalDegreesOfFreedom() > 0)
            {
                label idQ = 
                    cloud_.constProps(typeIdQ).charDissQuantumLevel()[0];
                label imaxQ = 0;
                
                //Only diatomics can be dissociated
                EcQ = translationalEnergy + EVibQ;

                imaxQ = EcQ/(physicoChemical::k.value()*cloud_.constProps(typeIdQ).thetaV()[0]);
                
                if(imaxQ-idQ > 0)
                {
                    //Dissociation can occur
                    totalReactionProbability += 1.0;
                    reactionProbabilities[0] = 1.0;
                }
            }
        }
        
        scalar ionisationEnergy = cloud_.constProps(typeIdQ).ionisationTemperature()*physicoChemical::k.value();
        
        // calculate if an ionisation of species Q is possible
        EcQ = translationalEnergy + EElecQ;

        if((EcQ - ionisationEnergy) > VSMALL)
        {
            //Ionisation can occur
            totalReactionProbability += 1.0;
            reactionProbabilities[1] = 1.0;
        }
        
        scalar heatOfReactionExchJoules = heatOfReactionExch_*physicoChemical::k.value();
        
        EcQ = translationalEnergy + EElecQ;
        
        scalar aDash = aCoeff_*(pow(2.5 - omegaPQ, bCoeff_)*exp(lgamma(2.5 - omegaPQ))/exp(lgamma(2.5 - omegaPQ + bCoeff_)));
         
        TColl = (translationalEnergy/(physicoChemical::k.value()))/(2.5 - omegaPQ);

        scalar activationEnergy = (aDash*pow((TColl/273.0) , bCoeff_) * fabs(heatOfReactionExchJoules));
        
//         scalar activationEnergy = (aCoeff_*pow((5000/273.0) , bCoeff_) * fabs(heatOfReactionExchJoules));

        if(heatOfReactionExchJoules < 0.0)
        {
            activationEnergy -= heatOfReactionExchJoules;
        }
        
        if(EcQ > activationEnergy)
        {
            label keyElectronicLevel = -1;
            
            for(label i = 0 ; i < maxElectronicLevelQ ; i++)
            {           
                scalar electronicEnergy = EElistQ[i];
                
                if(electronicEnergy > activationEnergy)
                {
                    break;
                } 
                
                keyElectronicLevel++;
            }
            
            EcQ = translationalEnergy + EElecQ /*+ heatOfReactionExchJoules*/;

            label trialELevel = cloud_.postCollisionElectronicEnergyLevel
                            (
                                EcQ,
                                maxElectronicLevelQ,
                                omegaPQ,
                                EElistQ,
                                gListQ
                            );
                            
            if(trialELevel == keyElectronicLevel)
            {
                scalar prob = 0.0;
                
                label nPossStates = 0;
                    
                if(maxElectronicLevelQ == 1)
                {
                    nPossStates = gListQ[0];
                }
                else
                {
                    forAll(EElistQ, n)
                    {
                        if(EcQ > EElistQ[n])
                        {
                            nPossStates += gListQ[n];
                        }
                    }
                }
                
                label nState = ceil(cloud_.rndGen().scalar01()*(nPossStates));
                label nAvailableStates = 0;
                label nLevel = -1;
                
                forAll(EElistQ, n)
                {
                    nAvailableStates += gListQ[n];
                    
                    if(nState <= nAvailableStates && nLevel < 0)
                    {
                        nLevel = n;
                    }
                }
                
                //Calculate the probability of it occuring
                
                scalar summation = 0.0;
                
                for(label i = 0 ; i <= nLevel; i++)
                { 
                    summation += gListQ[i]*pow((EcQ - EElistQ[i]),1.5-omegaPQ);
                }
                
                prob = (gListQ[trialELevel]*pow((EcQ - EElistQ[trialELevel]),1.5-omegaPQ))/summation;
                
                if(prob > cloud_.rndGen().scalar01())
                {
                    //Charge exchange can occur
                    totalReactionProbability += prob;
                    reactionProbabilities[2] = prob;
                }
            }
        }
        
        //Decide if a reaction is to occur
        
        if(totalReactionProbability > cloud_.rndGen().scalar01())
        {
            //A chemical reaction is to occur, choose which one
            
            scalarList normalisedProbabilities(reactionProbabilities.size(), 0.0);
            scalar cumulativeProbability = 0.0;
            
            normalisedProbabilities = reactionProbabilities/totalReactionProbability;
            
            forAll(normalisedProbabilities, i)
            {              
                //If current reaction can't occur, don't check for it
                if(normalisedProbabilities[i] > VSMALL)
                {
                    cumulativeProbability += normalisedProbabilities[i];
                    
                    if(cumulativeProbability > cloud_.rndGen().scalar01())
                    {
                        //Current reaction is to occur
                        
                        if(i == 0)
                        {
                            //Dissociation is to occur
                            dissocReactionQ = true;
                            break;
                        }
                        if(i == 1)
                        {
                            //Ionisation reaction is to occur
                            ionisationReactionQ = true;
                            break;
                        }
                        if(i == 2)
                        {
                            //Dissociation is to occur
                            chargeExchange = true;
                            break;
                        }
                    }
                }
            }
        }
        
        if(dissocReactionQ)
        {
            nDissociationReactions_++;
            nDissociationReactionsPerTimeStep_++;
            
            if(allowSplitting_)
            {
                relax_ = false;
                
                const scalar& heatOfReactionDissociationJoules = heatOfReactionDissoc_*physicoChemical::k.value();
                
                translationalEnergy = translationalEnergy + heatOfReactionDissociationJoules + EVibQ;
               
                scalar thetaVP = 0;
                scalar thetaDP = 0;
                scalar ZrefP = 0;
                scalar refTempZvP = 0;
                
                if(cloud_.constProps(typeIdP).vibrationalDegreesOfFreedom() > 0)
                {
                    thetaVP = cloud_.constProps(typeIdP).thetaV()[0];
                    thetaDP = cloud_.constProps(typeIdP).thetaD()[0];
                    ZrefP = cloud_.constProps(typeIdP).Zref()[0];
                    refTempZvP = cloud_.constProps(typeIdP).TrefZv()[0];
                }
                scalar jMaxP = cloud_.constProps(typeIdP).numberOfElectronicLevels();
                scalar rotationalDofP = cloud_.constProps(typeIdP).rotationalDegreesOfFreedom();

                translationalEnergy += EElecP;
                    
                label ELevelP = cloud_.postCollisionElectronicEnergyLevel
                                (
                                    translationalEnergy,
                                    jMaxP,
                                    omegaPQ,
                                    EElistP,
                                    gListP
                                );
                                
                translationalEnergy -= EElistP[ELevelP];
                
                label vibLevelP = 0;
                
                if(rotationalDofP > VSMALL)
                {
                    translationalEnergy += EVibP;
                    
                    label iMax = (translationalEnergy / (physicoChemical::k.value()*thetaVP));
                    
                    vibLevelP = cloud_.postCollisionVibrationalEnergyLevel
                                    (
                                            true,
                                            p.vibLevel()[0],
                                            iMax,
                                            thetaVP,
                                            thetaDP,
                                            refTempZvP,
                                            omegaPQ,
                                            ZrefP,
                                            translationalEnergy
                                        );
                                    
                    translationalEnergy -= vibLevelP*thetaVP*physicoChemical::k.value();
                                    
                    translationalEnergy += ERotP;
                    
                    ERotP = translationalEnergy*cloud_.postCollisionRotationalEnergy(rotationalDofP,ChiB);
                            
                    translationalEnergy -= ERotP;
                }
                
                scalar relVelNonDissoMol = sqrt(2.0*translationalEnergy/mR);

                //center of mass velocity of all particles
                vector Ucm = (mP*UP + mQ*UQ)/(mP + mQ);    

                // Variable Hard Sphere collision part

                scalar cosTheta = 2.0*cloud_.rndGen().scalar01() - 1.0;
            
                scalar sinTheta = sqrt(1.0 - cosTheta*cosTheta);
            
                scalar phi = twoPi*cloud_.rndGen().scalar01();
            
                vector postCollisionRelU =
                    relVelNonDissoMol
                    *vector
                        (
                            cosTheta,
                            sinTheta*cos(phi),
                            sinTheta*sin(phi)
                        );

                UP = Ucm + (postCollisionRelU*mQ/(mP + mQ)); 
                UQ = Ucm - (postCollisionRelU*mP/(mP + mQ));

                const label& typeId1 = dissociationProductIds_[0];
                const label& typeId2 = dissociationProductIds_[1];
                
                //Mass of Product one and two
                scalar mP1 = cloud_.constProps(typeId1).mass();
                scalar mP2 = cloud_.constProps(typeId2).mass();
                
                scalar mRatoms = mP1*mP2/(mP1 + mP2);
                
                translationalEnergy = ERotQ + EElecQ;
                
                scalar cRatoms = sqrt(2.0*translationalEnergy/mRatoms);

                // Variable Hard Sphere collision part
                scalar cosTheta2 = 2.0*cloud_.rndGen().scalar01() - 1.0;
            
                scalar sinTheta2 = sqrt(1.0 - cosTheta2*cosTheta2);
            
                scalar phi2 = twoPi*cloud_.rndGen().scalar01();
            
                vector postCollisionRelU2 = cRatoms
                *vector
                    (
                        cosTheta2,
                        sinTheta2*cos(phi2),
                        sinTheta2*sin(phi2)
                    );


                vector uQ1 = UQ + postCollisionRelU2*mP2/(mP1 + mP2);
                vector uQ2 = UQ - postCollisionRelU2*mP1/(mP1 + mP2);

                // P remains NON-DISSOCIATED.
                p.U() = UP;
                p.ERot() = ERotP;
                if(rotationalDofP > VSMALL)
                {
                    p.vibLevel().setSize(1, vibLevelP);
                }
                else 
                {
                    p.vibLevel().setSize(0, vibLevelP);
                }
                p.ELevel() = ELevelP;

                // Molecule Q will dissociate.
                vector position = q.position();
                
                label cell = -1;
                label tetFace = -1;
                label tetPt = -1;

                mesh_.findCellFacePt
                (
                    position,
                    cell,
                    tetFace,
                    tetPt
                );
                
                q.typeId() = typeId2;
                q.U() = uQ1;
                q.vibLevel().setSize(0,0);
                q.ERot() = 0.0;
                q.ELevel() = 0;
                
                label classificationQ = q.classification();
                scalar RWF = q.RWF();
                labelList vibLevel;
                
                // insert new product 2
                cloud_.addNewParcel
                (
                    position,
                    uQ2,
                    RWF,
                    0.0,
                    0,
                    cell,
                    tetFace,
                    tetPt,
                    typeId1,
                    -1,
                    classificationQ,
                    vibLevel
                );
            }
        }
        
        if(ionisationReactionQ)
        {
            nIonisationReactions_++;
            nIonisationReactionsPerTimeStep_++;
            
            if(allowSplitting_)
            {
                relax_ = false;
                
                const scalar& heatOfReactionIonisationJoules = heatOfReactionIon_*physicoChemical::k.value();
                
                translationalEnergy = translationalEnergy + heatOfReactionIonisationJoules + EElecQ;
                
                scalar thetaVP = 0;
                scalar thetaDP = 0;
                scalar ZrefP = 0;
                scalar refTempZvP = 0;
                
                if(cloud_.constProps(typeIdP).vibrationalDegreesOfFreedom() > 0)
                {
                    thetaVP = cloud_.constProps(typeIdP).thetaV()[0];
                    thetaDP = cloud_.constProps(typeIdP).thetaD()[0];
                    ZrefP = cloud_.constProps(typeIdP).Zref()[0];
                    refTempZvP = cloud_.constProps(typeIdP).TrefZv()[0];
                }
                
                scalar jMaxP = cloud_.constProps(typeIdP).numberOfElectronicLevels();
                scalar rotationalDofP = cloud_.constProps(typeIdP).rotationalDegreesOfFreedom();
                scalar rotationalDofQ = cloud_.constProps(typeIdQ).rotationalDegreesOfFreedom();

                translationalEnergy += EElecP;
                    
                label ELevelP = cloud_.postCollisionElectronicEnergyLevel
                                (
                                    translationalEnergy,
                                    jMaxP,
                                    omegaPQ,
                                    EElistP,
                                    gListP
                                );
                                
                translationalEnergy -= EElistP[ELevelP];
                
                label vibLevelP = 0;
                
                if(rotationalDofP > VSMALL)
                {
                    translationalEnergy += EVibP;
                    
                    label iMax = (translationalEnergy / (physicoChemical::k.value()*thetaVP));
                    
                    vibLevelP = cloud_.postCollisionVibrationalEnergyLevel
                                    (
                                            true,
                                            p.vibLevel()[0],
                                            iMax,
                                            thetaVP,
                                            thetaDP,
                                            refTempZvP,
                                            omegaPQ,
                                            ZrefP,
                                            translationalEnergy
                                        );
                                    
                    translationalEnergy -= vibLevelP*thetaVP*physicoChemical::k.value();
                                    
                    translationalEnergy += ERotP;
                    
                    ERotP = translationalEnergy*cloud_.postCollisionRotationalEnergy(rotationalDofP,ChiB);
                            
                    translationalEnergy -= ERotP;
                }
                
                scalar relVelNonDissoMol = sqrt(2.0*translationalEnergy/mR);

                //center of mass velocity of all particles
                vector Ucm = (mP*UP + mQ*UQ)/(mP + mQ);    

                // Variable Hard Sphere collision part

                scalar cosTheta = 2.0*cloud_.rndGen().scalar01() - 1.0;
            
                scalar sinTheta = sqrt(1.0 - cosTheta*cosTheta);
            
                scalar phi = twoPi*cloud_.rndGen().scalar01();
            
                vector postCollisionRelU =
                    relVelNonDissoMol
                    *vector
                        (
                            cosTheta,
                            sinTheta*cos(phi),
                            sinTheta*sin(phi)
                        );

                UP = Ucm + (postCollisionRelU*mQ/(mP + mQ)); 
                UQ = Ucm - (postCollisionRelU*mP/(mP + mQ));

                const label& typeId1 = ionisationProductIds_[0];
                const label& typeId2 = ionisationProductIds_[1];
                
                //Mass of Product one and two
                scalar mP1 = cloud_.constProps(typeId1).mass();
                scalar mP2 = cloud_.constProps(typeId2).mass();
                
                scalar mRatoms = mP1*mP2/(mP1 + mP2);
                
                translationalEnergy = ERotQ + EVibQ;
                
                scalar cRatoms = sqrt(2.0*translationalEnergy/mRatoms);

                // Variable Hard Sphere collision part
                scalar cosTheta2 = 2.0*cloud_.rndGen().scalar01() - 1.0;
            
                scalar sinTheta2 = sqrt(1.0 - cosTheta2*cosTheta2);
            
                scalar phi2 = twoPi*cloud_.rndGen().scalar01();
            
                vector postCollisionRelU2 = cRatoms
                *vector
                    (
                        cosTheta2,
                        sinTheta2*cos(phi2),
                        sinTheta2*sin(phi2)
                    );


                vector uQ1 = UQ + postCollisionRelU2*mP2/(mP1 + mP2);
                vector uQ2 = UQ - postCollisionRelU2*mP1/(mP1 + mP2);

                // P remains NON-DISSOCIATED.
                p.U() = UP;
                p.ERot() = ERotP;
                if(rotationalDofP > VSMALL)
                {
                    p.vibLevel().setSize(1, vibLevelP);
                }
                else
                {
                    p.vibLevel().setSize(0,0);
                }
                p.ELevel() = ELevelP;

                // Molecule Q will ionise.
                vector position = q.position();
                
                label cell = -1;
                label tetFace = -1;
                label tetPt = -1;

                mesh_.findCellFacePt
                (
                    position,
                    cell,
                    tetFace,
                    tetPt
                );
                
                q.typeId() = typeId1;
                q.U() = uQ1;
                if(rotationalDofQ > VSMALL)
                {
                    q.vibLevel().setSize(1,0);
                }       
                else
                {
                    q.vibLevel().setSize(0,0);
                } 
                q.ERot() = 0.0;
                q.ELevel() = 0;
                
                label classificationQ = q.classification();
                scalar RWF = q.RWF();
                labelList vibLevel;
                
                // insert new product 2
                cloud_.addNewParcel
                (
                    position,
                    uQ2,
                    RWF,
                    0.0,
                    0,
                    cell,
                    tetFace,
                    tetPt,
                    typeId2,
                    -1,
                    classificationQ,
                    vibLevel
                );
            }
        }
        
        if(chargeExchange)
        {
            nChargeExchangeReactions_++;
            nChargeExchangeReactionsPerTimeStep_++;
                    
            if(allowSplitting_)
            {
                relax_ = false;
                
                translationalEnergy = translationalEnergy + heatOfReactionExchJoules;
                
                translationalEnergy += ERotP + EVibP + EElecP + ERotQ + EVibQ + EElecQ;
                
                scalar relVel = sqrt((2.0*translationalEnergy)/mR);
                    
                    // centre of mass velocity of molecules (pre-split)
                vector Ucm = (mP*UP + mQ*UQ)/(mP + mQ);

                // Variable Hard Sphere collision part for collision of molecules
                scalar cosTheta = 2.0*cloud_.rndGen().scalar01() - 1.0;
            
                scalar sinTheta = sqrt(1.0 - cosTheta*cosTheta);
            
                scalar phi = twoPi*cloud_.rndGen().scalar01();
            
                vector postCollisionRelU =
                    relVel
                    *vector
                        (
                            cosTheta,
                            sinTheta*cos(phi),
                            sinTheta*sin(phi)
                        );
                        
                mP = cloud_.constProps(chargeExchangeProductIds_[0]).mass();
                mQ = cloud_.constProps(chargeExchangeProductIds_[1]).mass();
        
                UP = Ucm + (postCollisionRelU*mQ/(mP + mQ));
                UQ = Ucm - (postCollisionRelU*mP/(mP + mQ));
                
                p.typeId() = chargeExchangeProductIds_[1];
                p.U() = UP;
                p.ERot() = 0.0;
                 
                if(cloud_.constProps(chargeExchangeProductIds_[1]).
                                        rotationalDegreesOfFreedom() 
                    > VSMALL)
                {
                    p.vibLevel().setSize(1,0);
                }       
                else
                {
                    p.vibLevel().setSize(0,0);
                } 
                p.ELevel() = 0;
                
                q.typeId() = chargeExchangeProductIds_[0];
                q.U() = UQ;
                q.ERot() = 0.0;
                 
                if(cloud_.constProps(chargeExchangeProductIds_[0]).
                                        rotationalDegreesOfFreedom() 
                    > VSMALL)
                {
                    q.vibLevel().setSize(1,0);
                }       
                else
                {
                    q.vibLevel().setSize(0,0);
                } 
                q.ELevel() = 0;
            }
        }
    }
}

void  chargeExchange::outputResults(const label& counterIndex)
{
    if(writeRatesToTerminal_ == true)
    {
        // measure density 

        const List< DynamicList<dsmcParcel*> >& cellOccupancy
            = cloud_.cellOccupancy();
            
        volume_ = 0.0;

        forAll(cellOccupancy, c)
        {
            volume_ += mesh_.cellVolumes()[c];
        }

        List<label> mols;
	mols.append(0); mols.append(0);
        scalar volume = volume_;
        label nTotChargeExchangeReactions = nChargeExchangeReactions_;
        label nTotDissociationReactions = nDissociationReactions_;
        label nTotIonisationReactions = nIonisationReactions_;

        forAll(cellOccupancy, c)
        {
            const List<dsmcParcel*>& parcelsInCell = cellOccupancy[c];

            forAll(parcelsInCell, pIC)
            {
                dsmcParcel* p = parcelsInCell[pIC];

                label id = findIndex(reactantIds_, p->typeId());

                if(id != -1)
                {
                    mols[id]++;
                }
            }
        }

        //- Parallel communication
        if(Pstream::parRun())
        {
            reduce(mols[0], sumOp<label>());
            reduce(mols[1], sumOp<label>());
            reduce(volume, sumOp<scalar>());
            reduce(nTotChargeExchangeReactions, sumOp<label>());
            reduce(nTotDissociationReactions, sumOp<label>());
            reduce(nTotIonisationReactions, sumOp<label>());
        }

        numberDensities_[0] = (mols[0]*cloud().nParticle())/volume;
        
        if(reactantIds_[0] == reactantIds_[1])
        {
            numberDensities_[1] = (mols[0]*cloud().nParticle())/volume;
        }
        else
        {
            numberDensities_[1] = (mols[1]*cloud().nParticle())/volume;
        }

        const scalar& deltaT = mesh_.time().deltaT().value();

        word reactantMolA = cloud_.typeIdList()[reactantIds_[0]];
        word reactantMolB = cloud_.typeIdList()[reactantIds_[1]];

        word productMolA = cloud_.typeIdList()[chargeExchangeProductIds_[0]];
        word productMolB = cloud_.typeIdList()[chargeExchangeProductIds_[1]];
        
        word productMolC;
        word productMolD;
        
        if(dissociationPossible_)
        {
            productMolC = cloud_.typeIdList()[dissociationProductIds_[0]];
            productMolD = cloud_.typeIdList()[dissociationProductIds_[1]];
        }
        
        word productMolE = cloud_.typeIdList()[ionisationProductIds_[0]];
        word productMolF = cloud_.typeIdList()[ionisationProductIds_[1]];

        if((numberDensities_[0] > 0.0) && (numberDensities_[1] > 0.0))
        {
            scalar reactionRateChargeExchange = 0.0;
            scalar reactionRateDissociation = 0.0;
            scalar reactionRateIonisation = 0.0;
            
            reactionRateChargeExchange =
            (
                nTotChargeExchangeReactions
                * cloud_.nParticle()
            )/(counterIndex*deltaT*numberDensities_[0]* numberDensities_[1]*volume);
            
            Info<< "Charge exchange reaction " 
                <<  reactantMolA << " + " << reactantMolB
                <<  " --> " 
                << productMolA << " + " << productMolB
                << ", reaction rate = " << reactionRateChargeExchange
                << endl;
                
            reactionRateIonisation =
            (
                nTotIonisationReactions
                * cloud_.nParticle()
            )/(counterIndex*deltaT*numberDensities_[0]* numberDensities_[1]*volume);
            
            Info<< "Ionisation reaction " 
                <<  reactantMolA << " + " << reactantMolB
                <<  " --> " 
                << productMolE << " + " << productMolF
                << ", reaction rate = " << reactionRateIonisation
                << endl;
                
            if(dissociationPossible_)
            {
                reactionRateDissociation =
                (
                    nTotDissociationReactions
                    * cloud_.nParticle()
                )/(counterIndex*deltaT*numberDensities_[0]* numberDensities_[1]*volume);
                
                Info<< "Dissociation reaction " 
                    <<  reactantMolA << " + " << reactantMolB
                    <<  " --> " 
                    << productMolC << " + " << productMolD
                    << ", reaction rate = " << reactionRateDissociation
                    << endl;
            }
        }
    }
    else
    {
        label nTotChargeExchangeReactions = nChargeExchangeReactions_;  
        label nTotDissociationReactions = nDissociationReactions_;
        label nTotIonisationReactions = nIonisationReactions_;
        
        label nChargeExchangeReactionsPerTimeStep = nChargeExchangeReactionsPerTimeStep_;
        label nDissociationReactionsPerTimeStep = nDissociationReactionsPerTimeStep_;
        label nIonisationReactionsPerTimeStep = nIonisationReactionsPerTimeStep_;
        
        
        if(Pstream::parRun())
        {
            reduce(nTotChargeExchangeReactions, sumOp<label>());
            reduce(nTotDissociationReactions, sumOp<label>());
            reduce(nTotIonisationReactions, sumOp<label>());
            reduce(nChargeExchangeReactionsPerTimeStep, sumOp<label>());
            reduce(nDissociationReactionsPerTimeStep, sumOp<label>());
            reduce(nIonisationReactionsPerTimeStep, sumOp<label>());
        }
        
        if(nTotChargeExchangeReactions > VSMALL)
        {
                word reactantMolA = cloud_.typeIdList()[reactantIds_[0]];
                word reactantMolB = cloud_.typeIdList()[reactantIds_[1]];

                word productMolA = cloud_.typeIdList()[chargeExchangeProductIds_[0]];
                word productMolB = cloud_.typeIdList()[chargeExchangeProductIds_[1]];
            
                Info<< "Charge exchange reaction " 
                    <<  reactantMolA << " + " << reactantMolB
                    <<  " --> " 
                    << productMolA << " + " << productMolB
                    << " is active, nReactions this time step = " << nChargeExchangeReactionsPerTimeStep << endl;
        } 
        
        if(nTotDissociationReactions > VSMALL)
        {
                word reactantMolA = cloud_.typeIdList()[reactantIds_[0]];
                word reactantMolB = cloud_.typeIdList()[reactantIds_[1]];

                word productMolA = cloud_.typeIdList()[dissociationProductIds_[0]];
                word productMolB = cloud_.typeIdList()[dissociationProductIds_[1]];
            
                Info<< "Dissociation reaction " 
                    <<  reactantMolA << " + " << reactantMolB
                    <<  " --> " 
                    << productMolA << " + " << productMolB << " + " << reactantMolB
                    << " is active, nReactions this time step = " << nDissociationReactionsPerTimeStep << endl;
        } 
        
        if(nTotIonisationReactions > VSMALL)
        {
                word reactantMolA = cloud_.typeIdList()[reactantIds_[0]];
                word reactantMolB = cloud_.typeIdList()[reactantIds_[1]];

                word productMolA = cloud_.typeIdList()[ionisationProductIds_[0]];
                word productMolB = cloud_.typeIdList()[ionisationProductIds_[1]];
            
                Info<< "Ionisation reaction " 
                    <<  reactantMolA << " + " << reactantMolB
                    <<  " --> " 
                    << productMolA << " + " << productMolB << " + " << reactantMolB
                    << " is active, nReactions this time step = " << nIonisationReactionsPerTimeStep << endl;
        } 
    }

    nChargeExchangeReactionsPerTimeStep_ = 0.0;
    nDissociationReactionsPerTimeStep_ = 0.0;
    nIonisationReactionsPerTimeStep_ = 0.0;
}


const bool& chargeExchange::relax() const
{
    return relax_;
}

}
// End namespace Foam

// ************************************************************************* //
