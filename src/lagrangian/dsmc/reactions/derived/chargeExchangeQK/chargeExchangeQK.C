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

#include "chargeExchangeQK.H"
#include "addToRunTimeSelectionTable.H"
#include "fvc.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(chargeExchangeQK, 0);

addToRunTimeSelectionTable(dsmcReaction, chargeExchangeQK, dictionary);


// * * * * * * * * * * *  Protected Member functions * * * * * * * * * * * * //

void chargeExchangeQK::setProperties()
{
    dsmcReaction::setProperties();
    
    if (reactantIds_.size() != 2)
    {
        //- There must be exactly 2 reactants
        FatalErrorIn("chargeExchangeQK::setProperties()")
            << "For reaction named " << reactionName_ << nl
            << "There should be two reactants, instead of " 
            << reactantIds_.size() << nl 
            << exit(FatalError);
    }

    bool neutralFound = false;
    bool ionisedFound = false;
    
    forAll(reactantIds_, r)
    {
        switch (cloud_.constProps(reactantIds_[r]).charge())
        {
            case -1:
                //- No electron in the list of reactants
                FatalErrorIn("chargeExchangeQK::setProperties()")
                    << "For reaction named " << reactionName_ << nl
                    << "Reactant " << cloud_.typeIdList()[reactantIds_[r]]
                    << " is an electron and only a combination of " << nl 
                    << "one neutral and one ionised particle is accepted" << nl
                    << exit(FatalError);
                break;    
            case 0:
                neutralFound = true;
                break;
            case 1:
                ionisedFound = true;
                break;
        }
    }
        
    if (!neutralFound)
    {
        FatalErrorIn("chargeExchangeQK::setProperties()")
            << "For reaction named " << reactionName_ << nl
            << "None of the reactants is a neutral particle." << nl 
            << exit(FatalError);
    }
    else if (!ionisedFound)
    {
        FatalErrorIn("chargeExchangeQK::setProperties()")
            << "For reaction named " << reactionName_ << nl
            << "None of the reactants is an ionised particle." << nl 
            << exit(FatalError);
    }
    
    //- Reading in charge exchange products
    const wordList productsChargeExchange
    (
        propsDict_.lookup("chargeExchangeProducts")
    );

    if (productsChargeExchange.size() != 2)
    {
        FatalErrorIn("chargeExchangeQK::setProperties()")
            << "For reaction named " << reactionName_ << nl
            << "There should be two products, instead of " 
            << productsChargeExchange.size() << nl 
            << exit(FatalError);
    }
    
    productIdsChargeExchange_.setSize(productsChargeExchange.size());
    
    neutralFound = false;
    ionisedFound = false;

    forAll(productIdsChargeExchange_, r)
    {
        const label productIndex = 
            findIndex
            (
                cloud_.typeIdList(), 
                productsChargeExchange[r]
            );

        //- Check that products belong to the typeIdList as defined in 
        //  constant/dsmcProperties
        if (productIndex == -1)
        {
            FatalErrorIn("chargeExchangeQK::setProperties()")
                << "For reaction named " << reactionName_ << nl
                << "Cannot find type id: " << productsChargeExchange[r] << nl 
                << exit(FatalError);
        }

        switch (cloud_.constProps(productIndex).charge())
        {
            case -1:
                //- No electron in the list of products
                FatalErrorIn("chargeExchangeQK::setProperties()")
                    << "For reaction named " << reactionName_ << nl
                    << "Product " << cloud_.typeIdList()[reactantIds_[r]]
                    << " is an electron and only a combination of " << nl 
                    << "one neutral and one ionised particle is accepted" << nl
                    << exit(FatalError);
                break;    
            case 0:
                neutralFound = true;
                //- The neutral particle is set to be the first product
                productIdsChargeExchange_[0] = productIndex;
                break;
            case 1:
                ionisedFound = true;
                //- The ionised particle is set to be the second product
                productIdsChargeExchange_[1] = productIndex;
                break;
        }
    }
    
    if (!neutralFound)
    {
        FatalErrorIn("chargeExchangeQK::setProperties()")
            << "For reaction named " << reactionName_ << nl
            << "None of the products is a neutral particle." << nl 
            << exit(FatalError);
    }
    else if (!ionisedFound)
    {
        FatalErrorIn("chargeExchangeQK::setProperties()")
            << "For reaction named " << reactionName_ << nl
            << "None of the products is an ionised particle." << nl 
            << exit(FatalError);
    }
}


void chargeExchangeQK::testChargeExchange
(
    const dsmcParcel& p,
    const scalar translationalEnergy,
    const scalar omegaPQ,
    scalar& totalReactionProbability,
    scalar& reactionProbability
)
{
    const label typeIdP = p.typeId();
    
    //- Collision temperature: Eq.(10) of Bird's QK paper.
    const scalar TColl = translationalEnergy/physicoChemical::k.value()
        /(2.5 - omegaPQ); 
    
    const scalar aDash = 
        aCoeffChEx_
       *(
            pow(2.5 - omegaPQ, bCoeffChEx_)
           *exp(lgamma(2.5 - omegaPQ))
           /exp(lgamma(2.5 - omegaPQ + bCoeffChEx_))
        );

    scalar activationEnergy = 
        (
            aDash*pow(TColl/273.0, bCoeffChEx_)
           *fabs(heatOfReactionChargeExchangeJoules_)
        );
    
    if (heatOfReactionChargeExchangeJoules_ < 0.0) 
    {
        //- forward (endothermic) charge exchange reaction
        activationEnergy -= heatOfReactionChargeExchangeJoules_;
    }

    const label maxElectronicLevelP =
        cloud_.constProps(typeIdP).nElectronicLevels();
    const labelList& gListP = cloud_.constProps(typeIdP).electronicDegeneracyList();
    const scalarList& EElistP = cloud_.constProps(typeIdP).electronicEnergyList();
    
    const scalar EEleP = EElistP[p.ELevel()];
    
    //- Total collision energy
    const scalar EcP = translationalEnergy + EEleP;
    
    //- Condition for the charge exchange reaction to possibly occur
    if (EcP > activationEnergy)
    {
        label keyElectronicLevel = -1;
        
        for(label i=0; i<maxElectronicLevelP; i++)
        {           
            if (EElistP[i] > activationEnergy)
            {
                break;
            } 
            
            keyElectronicLevel++;
        }
        
        const label trialELevel = cloud_.postCollisionElectronicEnergyLevel
            (
                EcP,
                maxElectronicLevelP,
                omegaPQ,
                EElistP,
                gListP
            );
                        
        if (trialELevel == keyElectronicLevel)
        {
            scalar probChEx = 0.0;
            label nPossStates = 0;
                
            if (maxElectronicLevelP == 1)
            {
                nPossStates = gListP[0];
            }
            else
            {
                forAll(EElistP, n)
                {
                    if (EcP > EElistP[n])
                    {
                        nPossStates += gListP[n];
                    }
                }
            }
            
            label nState = ceil(cloud_.rndGen().sample01<scalar>()*(nPossStates));
            label nAvailableStates = 0;
            label nLevel = -1;
            
            forAll(EElistP, n)
            {
                nAvailableStates += gListP[n];
                
                if (nState <= nAvailableStates && nLevel < 0)
                {
                    nLevel = n;
                }
            }
            
            //- Calculate the probability of it occuring
            scalar summation = 0.0;
            
            for(label i=0; i<=nLevel; i++)
            { 
                summation += gListP[i]*pow(EcP - EElistP[i], 1.5-omegaPQ);
            }
            
            probChEx =
                (
                    gListP[trialELevel]
                   *pow(EcP - EElistP[trialELevel], 1.5-omegaPQ)
                )
               /summation;
            
            if (probChEx > cloud_.rndGen().sample01<scalar>())
            {
                //- Charge exchange can occur
                reactionProbability = probChEx;
                totalReactionProbability += reactionProbability;
            }
        }
    }
}


void chargeExchangeQK::chargeExchange
(
    dsmcParcel& p,
    dsmcParcel& q
)
{
    const label typeIdP = p.typeId();
    const label typeIdQ = q.typeId();
    
    nTotChargeExchangeReactions_++;
    nChargeExchangeReactionsPerTimeStep_++;
    
    if (allowSplitting_)
    {
        relax_ = false;
        
        vector UP = p.U();
        vector UQ = q.U();
        
        const scalar mP = cloud_.constProps(typeIdP).mass();
        const scalar mQ = cloud_.constProps(typeIdQ).mass();
        const scalar mR = mP*mQ/(mP + mQ);
        const scalar cRsqr = magSqr(UP - UQ);
        
        scalar translationalEnergy = 0.5*mR*cRsqr;
        
        //- Center of mass velocity (pre-charge exchange)
        const vector Ucm = (mP*UP + mQ*UQ)/(mP + mQ);

        const label typeIdNeutral = productIdsChargeExchange_[0];
        const label typeIdIonised = productIdsChargeExchange_[1];

        //- Change species properties
        const scalar mPChEx = cloud_.constProps(typeIdIonised).mass();
        const scalar mQChEx = cloud_.constProps(typeIdNeutral).mass();
        const scalar mRChEx = mPChEx*mQChEx/(mPChEx + mQChEx);
        
        const scalar EVibP = cloud_.constProps(typeIdP).eVib_tot(p.vibLevel());
        const scalar EEleP = cloud_.constProps(typeIdP).electronicEnergyList()[p.ELevel()];
        const scalar EVibQ = cloud_.constProps(typeIdQ).eVib_tot(q.vibLevel());
        const scalar EEleQ = cloud_.constProps(typeIdQ).electronicEnergyList()[q.ELevel()];
        
        //- Assumption: no energy redistribution for both particles
        //  All the energy is stored in the translational mode
        translationalEnergy += p.ERot() + EVibP + EEleP + q.ERot() + EVibQ
             + EEleQ + heatOfReactionChargeExchangeJoules_;
            
        const scalar relVelChExMol = sqrt(2.0*translationalEnergy/mRChEx);
        
        //- Variable Hard Sphere collision part for collision of molecules
        const scalar cosTheta = 2.0*cloud_.rndGen().sample01<scalar>() - 1.0;
        const scalar sinTheta = sqrt(1.0 - cosTheta*cosTheta);
        const scalar phi = twoPi*cloud_.rndGen().sample01<scalar>();
    
        const vector postCollisionRelU =
            relVelChExMol
           *vector
            (
                cosTheta,
                sinTheta*cos(phi),
                sinTheta*sin(phi)
            );
        
        UP = Ucm + postCollisionRelU*mQChEx/(mPChEx + mQChEx);
        UQ = Ucm - postCollisionRelU*mPChEx/(mPChEx + mQChEx);
        
        //- p is originally the neutral particle and becomes the ionised particle
        p.typeId() = typeIdIonised;
        p.U() = UP;
        p.ERot() = 0.0; 
        p.vibLevel().setSize
        (
            cloud_.constProps
            (
                typeIdIonised
            ).nVibrationalModes(),
            0
        );
        p.ELevel() = 0;
        
        //- q is originally the ionised particle and becomes the neutral particle
        q.typeId() = typeIdNeutral;
        q.U() = UQ;
        q.ERot() = 0.0;
        q.vibLevel().setSize
        (
            cloud_.constProps
            (
                typeIdNeutral
            ).nVibrationalModes(),
            0
        );
        q.ELevel() = 0;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
chargeExchangeQK::chargeExchangeQK
(
    Time& t,
    dsmcCloud& cloud,
    const dictionary& dict
)
:
    dsmcReaction(t, cloud, dict),
    propsDict_(dict.subDict(typeName + "Properties")),
    posMolReactant_(-1),
    productIdsChargeExchange_(),
    chargeExchangeStr_(word::null),
    nTotChargeExchangeReactions_(0),
    nChargeExchangeReactionsPerTimeStep_(0),
    heatOfReactionChargeExchangeJoules_
    (
        readScalar(propsDict_.lookup("heatOfReactionChargeExchange"))
       *physicoChemical::k.value()
    ),
    aCoeffChEx_(readScalar(propsDict_.lookup("aCoeff_ChEx"))),
    bCoeffChEx_(readScalar(propsDict_.lookup("bCoeff_ChEx"))),
    volume_(0.0)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

chargeExchangeQK::~chargeExchangeQK()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void chargeExchangeQK::initialConfiguration()
{
    setProperties();
    
    const word& reactantA = cloud_.typeIdList()[reactantIds_[0]];
    const word& reactantB = cloud_.typeIdList()[reactantIds_[1]];
    
    const word& productA = cloud_.typeIdList()[productIdsChargeExchange_[0]];
    const word& productB = cloud_.typeIdList()[productIdsChargeExchange_[1]];
    
    chargeExchangeStr_ = "Charge exchange reaction " + reactantA + " + " 
        + reactantB + " --> " + productA + " + " + productB;
}


bool chargeExchangeQK::tryReactMolecules
(
    const label& typeIdP,
    const label& typeIdQ
)
{
    //- Function used when setting the pair addressing matrix
    const label reactantPId = findIndex(reactantIds_, typeIdP);
    const label reactantQId = findIndex(reactantIds_, typeIdQ);

    //- If both indices were found in the list of reactants, there Ids will be
    //  different from -1
    if ((reactantPId != -1) && (reactantQId != -1))
    {
        //- Case of dissimilar species colliding (by definition)
        if (reactantPId != reactantQId)
        {
            return true;
        }
    }
        
    return false;
}


void chargeExchangeQK::reaction
(
    dsmcParcel& p,
    dsmcParcel& q,
    const DynamicList<label>& candidateList,
    const List<DynamicList<label> >& candidateSubList,
    const label& candidateP,
    const List<label>& whichSubCell
)
{}


void chargeExchangeQK::reaction(dsmcParcel& p, dsmcParcel& q)
{
    //- Reset the relax switch
    relax_ = true;
    
    const label typeIdP = p.typeId();
    const label typeIdQ = q.typeId();
    
    //- Charge exchange reaction P + Q+ --> P+ + Q
    //  If P is the first (neutral) reactant
    //  NB: Q is necessarily second otherwise this class would not have been
    //  selected
    if (cloud_.constProps(typeIdP).charge() == 0) 
    { 
        const scalar mP = cloud_.constProps(typeIdP).mass();
        const scalar mQ = cloud_.constProps(typeIdQ).mass();
        const scalar mR = mP*mQ/(mP + mQ);
        
        const scalar omegaPQ =
            0.5
            *(
                  cloud_.constProps(typeIdP).omega()
                + cloud_.constProps(typeIdQ).omega()
            );
        
        const scalar cRsqr = magSqr(p.U() - q.U());
        const scalar translationalEnergy = 0.5*mR*cRsqr;
        
        //- Possible reactions:
        // 1. Charge exchange reaction
        
        scalar totalReactionProbability = 0.0;
        scalarList reactionProbabilities(1, 0.0);
        
        testChargeExchange
        (
            p,
            translationalEnergy,
            omegaPQ,
            totalReactionProbability,
            reactionProbabilities[0]
        );
        
        //- Decide if a charge exchange reaction is to occur
        if (totalReactionProbability > cloud_.rndGen().sample01<scalar>())
        {
            chargeExchange(p, q);
        }
    }
    else
    {
        //- If P is the second reactant M, then switch arguments in this
        //  function and P will be first
        chargeExchangeQK::reaction(q, p);
    }
}

void chargeExchangeQK::outputResults(const label& counterIndex)
{
    if (writeRatesToTerminal_)
    {
        //- measure density 
        const List<DynamicList<dsmcParcel*>>& cellOccupancy = cloud_.cellOccupancy();
            
        volume_ = 0.0;

        labelList molsReactants(2, 0);

        forAll(cellOccupancy, c)
        {
            const List<dsmcParcel*>& parcelsInCell = cellOccupancy[c];

            forAll(parcelsInCell, pIC)
            {
                dsmcParcel* p = parcelsInCell[pIC];
                
                const label pos = findIndex(reactantIds_, p->typeId());

                if (pos != -1)
                {
                    molsReactants[pos]++;
                }
            }

            volume_ += mesh_.cellVolumes()[c];
        }
        
        scalar volume = volume_;
        if (Pstream::parRun())
        {
            reduce(volume, sumOp<scalar>());
        }
        
        scalarList numberDensities(2, cloud_.nParticle()/volume);
        numberDensities[0] *= molsReactants[0];
        numberDensities[1] *= molsReactants[1];
        
        label nTotChargeExchangeReactions = nTotChargeExchangeReactions_;
        label nChargeExchangeReactionsPerTimeStep = nChargeExchangeReactionsPerTimeStep_;

        const scalar deltaT = mesh_.time().deltaT().value();
        scalar factor = 0.0;
        
        if (reactantIds_[0] == reactantIds_[1] && numberDensities[0] > 0.0)
        {
            factor = cloud_.nParticle()/
                (
                    counterIndex*deltaT
                   *numberDensities[0]*numberDensities[0]
                   *volume
                );
        }
        else if (numberDensities[0] > 0.0 && numberDensities[1] > 0.0)
        {
            factor = cloud_.nParticle()/
                (
                    counterIndex*deltaT
                   *numberDensities[0]*numberDensities[1]
                   *volume
                );
        }
        
        if (Pstream::parRun())
        {
            //- Parallel communication
            reduce(molsReactants[0], sumOp<label>());
            reduce(molsReactants[1], sumOp<label>());
            reduce(nTotChargeExchangeReactions, sumOp<label>());
            reduce(nChargeExchangeReactionsPerTimeStep, sumOp<label>());
        }
        
        const scalar reactionRateChargeExchange = factor*nTotChargeExchangeReactions;
        
        Info<< chargeExchangeStr_
            << ", reaction rate = " << reactionRateChargeExchange
            << ", nReactions = " << nChargeExchangeReactionsPerTimeStep
            << endl;
    }
    else
    {
        label nTotChargeExchangeReactions = nTotChargeExchangeReactions_;   
        label nChargeExchangeReactionsPerTimeStep = nChargeExchangeReactionsPerTimeStep_;
        
        if (Pstream::parRun())
        {
            //- Parallel communication
            reduce(nTotChargeExchangeReactions, sumOp<label>());
            reduce(nChargeExchangeReactionsPerTimeStep, sumOp<label>());
        }
        
        if (nTotChargeExchangeReactions > 0)
        {
            Info<< chargeExchangeStr_
                << " is active, nReactions this time step = " 
                << nChargeExchangeReactionsPerTimeStep
                << endl;
         }
    }

    nChargeExchangeReactionsPerTimeStep_ = 0;
}

}
// End namespace Foam

// ************************************************************************* //
