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

#include "dissociationIonisationExchange.H"
#include "addToRunTimeSelectionTable.H"
#include "fvc.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(dissociationIonisationExchange, 0);

addToRunTimeSelectionTable(dsmcReaction, dissociationIonisationExchange, dictionary);



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
dissociationIonisationExchange::dissociationIonisationExchange
(
    Time& t,
    dsmcCloud& cloud,
    const dictionary& dict
)
:
    dsmcReaction(t, cloud, dict),
    propsDict_(dict.subDict(typeName + "Properties")),
    reactantIds_(),
    exchangeProductIds_(),
    chargeExchangeProductIds_(),
    dissociationProductIds_(),
    ionisationProductsIdsP_(),
    ionisationProductsIdsQ_(),
    chargedAtom_(Switch(propsDict_.lookup("chargedAtom"))),
    chargedMolecule_(Switch(propsDict_.lookup("chargedMolecule"))),
    chargeExchange_(Switch(propsDict_.lookup("chargeExchange"))),
    heatOfReactionDiss_(),
    heatOfReactionExch_(readScalar(propsDict_.lookup("heatOfReactionExch"))),
    heatOfReactionChargeExchange_(),
    heatOfReactionIonP_(),
    heatOfReactionIonQ_(),
    aCoeff_(readScalar(propsDict_.lookup("aCoeff"))),
    bCoeff_(readScalar(propsDict_.lookup("bCoeff"))),
    aCoeffCharge_(),
    bCoeffCharge_(),
    nTotExchangeReactions_(0),
    nTotChargeExchangeReactions_(),
    nTotDissociationReactions_(0),
    nTotIonisationReactionsP_(0),
    nTotIonisationReactionsQ_(0),
    nExchangeReactionsPerTimeStep_(0),
    nChargeExchangeReactionsPerTimeStep_(0),
    nDissociationReactionsPerTimeStep_(0),
    nIonisationReactionsPPerTimeStep_(0),
    nIonisationReactionsQPerTimeStep_(0),
    reactionName_(propsDict_.lookup("reactionName")),
    relax_(true),
    allowSplitting_(false),
    writeRatesToTerminal_(false),
    volume_(0.0),
    numberDensities_(2, 0.0)
{
    if(!chargedMolecule_)
    {
        heatOfReactionDiss_ = readScalar(propsDict_.lookup("heatOfReactionDiss"));
        heatOfReactionIonP_ = readScalar(propsDict_.lookup("heatOfReactionIonP"));
    }
    
    if(!chargedAtom_)
    {
        heatOfReactionIonQ_ = readScalar(propsDict_.lookup("heatOfReactionIonQ"));
    }
    
    if(chargeExchange_)
    {
        heatOfReactionChargeExchange_ = readScalar(propsDict_.lookup("heatOfReactionChargeExch"));
        aCoeffCharge_ = readScalar(propsDict_.lookup("aCoeffChargeExchange"));
        bCoeffCharge_ = readScalar(propsDict_.lookup("bCoeffChargeExchange"));
    }

}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

dissociationIonisationExchange::~dissociationIonisationExchange()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void dissociationIonisationExchange::initialConfiguration()
{
    setProperties();
}

void dissociationIonisationExchange::setProperties()
{
    // reading in reactants


    const List<word> reactantMolecules (propsDict_.lookup("reactantMolecules"));

    if(reactantMolecules.size() != 2)
    {
        FatalErrorIn("dissociationIonisationExchange::setProperties()")
            << "There should be two or more reactants, instead of " 
            << reactantMolecules.size() << nl 
            << exit(FatalError);
    }
    
    if(reactantMolecules[0] == reactantMolecules[1])
    {
        FatalErrorIn("dissociationIonisationExchange::setProperties()")
            << "Reactant molecules cannot be same species." << nl
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
            FatalErrorIn("dissociationIonisationExchange::setProperties()")
                << "Cannot find type id: " << reactantMolecules[r] << nl 
                << exit(FatalError);
        }
    }
    
    // check that the first reactant is a 'MOLECULE' 

    const label& rDof1 = cloud_.constProps(reactantIds_[0]).rotationalDegreesOfFreedom();

    if(rDof1 < 1)
    {
        FatalErrorIn("dissociationIonisationExchange::setProperties()")
            << "Reactant 1 must be a molecule (not an atom): " << reactantMolecules[0] 
            << nl 
            << exit(FatalError);
    }
    
    const label& vDof1 = cloud_.constProps(reactantIds_[0]).vibrationalDegreesOfFreedom();

    if(vDof1 > 1)
    {
         FatalErrorIn("dissociationIonisationExchange::setProperties()")
            << "Reactions are currently only implemented for monatomic and diatomic species"
            << " This is a polyatomic:" << reactantMolecules[0] 
            << nl 
            << exit(FatalError);
    }
    
    // check that the second reactant is an 'ATOM' 

    const label& rDof2 = cloud_.constProps(reactantIds_[1]).rotationalDegreesOfFreedom();

    if(rDof2 > 0)
    {
        FatalErrorIn("dissociationIonisationExchange::setProperties()")
            << "Reactant 2 must be an atom (not a molecule): " << reactantMolecules[1] 
            << nl 
            << exit(FatalError);
    }

    // reading in products

    const List<word> exchangeProductMolecules (propsDict_.lookup("productsOfExchangeReaction"));

    if(exchangeProductMolecules.size() != 2)
    {
        FatalErrorIn("dissociationIonisationExchange::setProperties()")
            << "There should be two dissociationIonisationExchange reaction products, instead of " 
            << reactantMolecules.size() << nl 
            << exit(FatalError);
    }
    
    if(exchangeProductMolecules[0] == exchangeProductMolecules[1])
    {
        FatalErrorIn("dissociationIonisationExchange::setProperties()")
            << "Exchange reaction product molecules cannot be same species." << nl
            << exit(FatalError);
    }

    exchangeProductIds_.setSize(exchangeProductMolecules.size(), -1);

    forAll(exchangeProductIds_, r)
    {
        exchangeProductIds_[r] = findIndex(cloud_.typeIdList(), exchangeProductMolecules[r]);

        // check that reactants belong to the typeIdList (constant/dsmcProperties)
        if(exchangeProductIds_[r] == -1)
        {
            FatalErrorIn("dissociationIonisationExchange::setProperties()")
                << "Cannot find type id: " << exchangeProductMolecules[r] << nl 
                << exit(FatalError);
        }
    }
    
    //check that first exchange product is a 'MOLECULE' (not an 'ATOM')

    const scalar& rDof3 = cloud_.constProps(exchangeProductIds_[0]).rotationalDegreesOfFreedom();

    if(rDof3 < 0)
    {
        FatalErrorIn("dissociationIonisationExchange::setProperties()")
            << "First product of the exchange reaction must be a molecule (not an atom): " << exchangeProductMolecules[0] 
            << nl 
            << exit(FatalError);
    }

    //check that second exchange product is an 'ATOM' (not a 'MOLECULE')

    const scalar& rDof4 = cloud_.constProps(exchangeProductIds_[1]).rotationalDegreesOfFreedom();

    if(rDof4 > 0)
    {
        FatalErrorIn("dissociationIonisationExchange::setProperties()")
            << "Second product of the exchange reaction must be an atom (not a molecule): " << exchangeProductMolecules[1] 
            << nl 
            << exit(FatalError);
    }
    
    if(chargeExchange_)
    {
        // reading in products

        const List<word> chargeExchangeProductMolecules (propsDict_.lookup("productsOfChargeExchangeReaction"));

        if(chargeExchangeProductMolecules.size() != 2)
        {
            FatalErrorIn("dissociationIonisationExchange::setProperties()")
                << "There should be two charge exchange reaction products, instead of " 
                << reactantMolecules.size() << nl 
                << exit(FatalError);
        }
        
        if(chargeExchangeProductMolecules[0] == chargeExchangeProductMolecules[1])
        {
            FatalErrorIn("dissociationIonisationExchange::setProperties()")
                << "Charge exchange reaction product molecules cannot be same species." << nl
                << exit(FatalError);
        }

        chargeExchangeProductIds_.setSize(chargeExchangeProductMolecules.size(), -1);

        forAll(chargeExchangeProductIds_, r)
        {
            chargeExchangeProductIds_[r] = findIndex(cloud_.typeIdList(), chargeExchangeProductMolecules[r]);

            // check that reactants belong to the typeIdList (constant/dsmcProperties)
            if(chargeExchangeProductIds_[r] == -1)
            {
                FatalErrorIn("dissociationIonisationExchange::setProperties()")
                    << "Cannot find type id: " << chargeExchangeProductMolecules[r] << nl 
                    << exit(FatalError);
            }
        }
        
        //check that first exchange product is a 'MOLECULE' (not an 'ATOM')

        const scalar& rDof9 = cloud_.constProps(chargeExchangeProductIds_[0]).rotationalDegreesOfFreedom();

        if(rDof9 < 0)
        {
            FatalErrorIn("dissociationIonisationExchange::setProperties()")
                << "First product of the exchange reaction must be a molecule (not an atom): " << chargeExchangeProductMolecules[0] 
                << nl 
                << exit(FatalError);
        }

        //check that second exchange product is an 'ATOM' (not a 'MOLECULE')

        const scalar& rDof10 = cloud_.constProps(chargeExchangeProductIds_[1]).rotationalDegreesOfFreedom();

        if(rDof10 > 0)
        {
            FatalErrorIn("dissociationIonisationExchange::setProperties()")
                << "Second product of the exchange reaction must be an atom (not a molecule): " << chargeExchangeProductMolecules[1] 
                << nl 
                << exit(FatalError);
        }        
    }
    
    if(!chargedMolecule_)
    {
        const List<word> dissociationProductMolecules (propsDict_.lookup("productsOfDissociatedMolecule"));
        
        if(dissociationProductMolecules.size() != 2)
        {
            FatalErrorIn("dissociationIonisationExchange::setProperties()")
                << "There should be two dissociation reaction products, instead of " 
                << reactantMolecules.size() << nl 
                << exit(FatalError);
        }

        dissociationProductIds_.setSize(dissociationProductMolecules.size(), -1);

        forAll(dissociationProductIds_, r)
        {
            dissociationProductIds_[r] = findIndex(cloud_.typeIdList(), dissociationProductMolecules[r]);

            // check that reactants belong to the typeIdList (constant/dsmcProperties)
            if(dissociationProductIds_[r] == -1)
            {
                FatalErrorIn("dissociationIonisationExchange::setProperties()")
                    << "Cannot find type id: " << dissociationProductMolecules[r] << nl 
                    << exit(FatalError);
            }
        }
        
        const scalar& rDof5= cloud_.constProps(dissociationProductIds_[0]).rotationalDegreesOfFreedom();

        if(rDof5 > 0)
        {
            FatalErrorIn("dissociationIonisationExchange::setProperties()")
                << "First product of the dissociation reaction must be an atom (not a molecule): " << dissociationProductMolecules[0] 
                << nl 
                << exit(FatalError);
        }
        
        //check that second exchange product is an 'ATOM' (not a 'MOLECULE')
        
        const scalar& rDof6 = cloud_.constProps(dissociationProductIds_[1]).rotationalDegreesOfFreedom();

        if(rDof6 > 0)
        {
            FatalErrorIn("dissociationIonisationExchange::setProperties()")
                << "Second product of the exchange reaction must be an atom (not a molecule): " << dissociationProductMolecules[1] 
                << nl 
                << exit(FatalError);
        }
        
        //read in ionisation products
        
        const List<word> ionisationProductMoleculesP (propsDict_.lookup("productsOfIonisedMolecule"));
        
        if(ionisationProductMoleculesP.size() != 2)
        {
            FatalErrorIn("dissociationIonisationExchange::setProperties()")
                << "There should be two ionisation reaction products, instead of " 
                << reactantMolecules.size() << nl 
                << exit(FatalError);
        }

        ionisationProductsIdsP_.setSize(ionisationProductMoleculesP.size(), -1);

        forAll(ionisationProductsIdsP_, r)
        {
            ionisationProductsIdsP_[r] = findIndex(cloud_.typeIdList(), ionisationProductMoleculesP[r]);

            // check that reactants belong to the typeIdList (constant/dsmcProperties)
            if(ionisationProductsIdsP_[r] == -1)
            {
                FatalErrorIn("dissociationIonisationExchange::setProperties()")
                    << "Cannot find type id: " << ionisationProductMoleculesP[r] << nl 
                    << exit(FatalError);
            }
        }
        
        //check that first ionisation product is a 'MOLECULE' (not an 'ATOM')
        
        const scalar& rDof7 = cloud_.constProps(ionisationProductsIdsP_[0]).rotationalDegreesOfFreedom();

        if(rDof7 < 0)
        {
            FatalErrorIn("dissociationIonisationExchange::setProperties()")
                << "First product of the molecule ionisation reaction must be a molecule (not an atom): " << ionisationProductMoleculesP[0] 
                << nl 
                << exit(FatalError);
        }
        
        //check that second ionisation product is a 'ELECTRON'
        
        const label& charge = cloud_.constProps(ionisationProductsIdsP_[1]).charge();

        if(charge != -1)
        {
            FatalErrorIn("dissociationIonisationExchange::setProperties()")
                << "Second product of the molecule ionisation reaction must be an electron: " << ionisationProductMoleculesP[1] 
                << nl 
                << exit(FatalError);
        }
    }
    
    if(!chargedAtom_)
    {
        const List<word> ionisationProductMoleculesQ (propsDict_.lookup("productsOfIonisedAtom"));
        
        if(ionisationProductMoleculesQ.size() != 2)
        {
            FatalErrorIn("dissociationIonisationExchange::setProperties()")
                << "There should be two ionisation reaction products, instead of " 
                << reactantMolecules.size() << nl 
                << exit(FatalError);
        }

        ionisationProductsIdsQ_.setSize(ionisationProductMoleculesQ.size(), -1);

        forAll(ionisationProductsIdsQ_, r)
        {
            ionisationProductsIdsQ_[r] = findIndex(cloud_.typeIdList(), ionisationProductMoleculesQ[r]);

            // check that reactants belong to the typeIdList (constant/dsmcProperties)
            if(ionisationProductsIdsQ_[r] == -1)
            {
                FatalErrorIn("dissociationIonisationExchange::setProperties()")
                    << "Cannot find type id: " << ionisationProductMoleculesQ[r] << nl 
                    << exit(FatalError);
            }
        }
        
        //check that first ionisation product is an 'ATOM'
        
        const scalar& rDof8 = cloud_.constProps(ionisationProductsIdsQ_[0]).rotationalDegreesOfFreedom();

        if(rDof8 > 0)
        {
            FatalErrorIn("dissociationIonisationExchange::setProperties()")
                << "First product of the atom ionisation reaction must be an atom (not a molecule): " << ionisationProductMoleculesQ[0] 
                << nl 
                << exit(FatalError);
        }
        
        //check that second ionisation product is a 'ELECTRON'
        
        const label& charge = cloud_.constProps(ionisationProductsIdsQ_[1]).charge();

        if(charge != -1)
        {
            FatalErrorIn("dissociationIonisationExchange::setProperties()")
                << "Second product of the atom ionisation reaction must be an electron: " << ionisationProductMoleculesQ[1] 
                << nl 
                << exit(FatalError);
        }
    }
}

bool dissociationIonisationExchange::tryReactMolecules(const label& typeIdP, const label& typeIdQ)
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

void dissociationIonisationExchange::reaction
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


void dissociationIonisationExchange::reaction
(
    dsmcParcel& p,
    dsmcParcel& q
)
{

    label typeIdP = p.typeId();
    label typeIdQ = q.typeId();

    //if particle p is the molecule and q is the atom...
    if(typeIdP == reactantIds_[0] && typeIdQ == reactantIds_[1])
    {       
        relax_ = true;
        
        scalar totalReactionProbability = 0.0;
        scalarList reactionProbabilities(5, 0.0);
        
        vector UP = p.U();
        vector UQ = q.U();
        scalar ERotP = p.ERot();
        scalar ERotQ = q.ERot();

        scalar EVibP = p.vibLevel()[0]*cloud_.constProps(typeIdP).thetaV()[0]*physicoChemical::k.value();

        scalar EEleP = cloud_.constProps(typeIdP).electronicEnergyList()[p.ELevel()];
        scalar EEleQ = cloud_.constProps(typeIdQ).electronicEnergyList()[q.ELevel()];

        scalar mP = cloud_.constProps(typeIdP).mass();
        scalar mQ = cloud_.constProps(typeIdQ).mass();

        scalar mR = mP*mQ/(mP + mQ);
        scalar cRsqr = magSqr(UP - UQ);
        scalar translationalEnergy = 0.5*mR*cRsqr;
        
        scalar thetaVP = cloud_.constProps(typeIdP).thetaV()[0];
        scalar thetaDP = cloud_.constProps(typeIdP).thetaD()[0];
        scalar ZrefP = cloud_.constProps(typeIdP).Zref()[0];
        scalar refTempZvP = cloud_.constProps(typeIdP).TrefZv()[0];
       
        label jMaxP = cloud_.constProps(typeIdP).numberOfElectronicLevels();
        List<label> gListP = cloud_.constProps(typeIdP).degeneracyList();
        List<scalar> EElistP = cloud_.constProps(typeIdP).electronicEnergyList();
       
        label jMaxQ = cloud_.constProps(typeIdQ).numberOfElectronicLevels();
        List<label> gListQ = cloud_.constProps(typeIdQ).degeneracyList();
        List<scalar> EElistQ = cloud_.constProps(typeIdQ).electronicEnergyList();
       
        scalar rotationalDofP = cloud_.constProps(typeIdP).rotationalDegreesOfFreedom();

        scalar omegaPQ =
        0.5
        *(
            cloud_.constProps(typeIdP).omega()
            + cloud_.constProps(typeIdQ).omega()
        );
        
        scalar ChiB = 2.5 - omegaPQ;

        scalar EcPQ = 0.0;
        scalar TColl = 0.0;
        label idP = -1;
        label deltaDissoIP= 0;
        label imaxP = 0;
        label iaP = 0;
        bool dissocReaction = false;
        bool ionisationReactionP = false;
        bool ionisationReactionQ = false;
        bool exchangeReaction = false;
        bool chargeExchangeReaction = false;
       
        if(!chargedMolecule_)
        {
            // firstly calculate dissociation probability (0 or 1).
            EcPQ = translationalEnergy + EVibP;

            imaxP = EcPQ/(physicoChemical::k.value()*cloud_.constProps(typeIdP).thetaV()[0]);

            idP = cloud_.constProps(typeIdP).thetaD()[0]/cloud_.constProps(typeIdP).thetaV()[0];

            deltaDissoIP = imaxP-idP;

            if(deltaDissoIP > 0)
            {
                // DISSOCIATION CAN OCCUR //
                
                totalReactionProbability += 1.0;
                reactionProbabilities[0] = 1.0;
            }
            
            //Now, ionisation of the molecule (P)
            
            scalar ionisationEnergy = cloud_.constProps(typeIdP).ionisationTemperature()*physicoChemical::k.value();
            
            // calculate if an ionisation of species P is possible
            EcPQ = translationalEnergy + EEleP;

            if((EcPQ - ionisationEnergy) > VSMALL)
            {
                //Ionisation can occur
                totalReactionProbability += 1.0;
                reactionProbabilities[1] = 1.0;
            }
        }
        
        
        
        if(!chargedAtom_)
        {
            //Now, ionisation of the atom (Q)
            
            scalar ionisationEnergy = cloud_.constProps(typeIdQ).ionisationTemperature()*physicoChemical::k.value();
            
            // calculate if an ionisation of species Q is possible
            EcPQ = translationalEnergy + EEleQ;

            if((EcPQ - ionisationEnergy) > VSMALL)
            {
                //Ionisation can occur
                totalReactionProbability += 1.0;
                reactionProbabilities[2] = 1.0;
            }
        }
        
        
        
        EcPQ = translationalEnergy + EVibP;
        
        scalar P_exch = 0.0;

        // Next, calculate exchange probability.

        TColl = (translationalEnergy/(physicoChemical::k.value()))/(2.5 - omegaPQ);// Eq.(10) of Bird QK paper.
            
        scalar heatOfReactionExchJoules = heatOfReactionExch_*physicoChemical::k.value();
        
        scalar aDash = aCoeff_*(pow(2.5 - omegaPQ, bCoeff_)*exp(lgamma(2.5 - omegaPQ))/exp(lgamma(2.5 - omegaPQ + bCoeff_)));

        scalar activationEnergy = (aDash*pow((TColl/273.0) , bCoeff_) * fabs(heatOfReactionExchJoules));
        
        if(heatOfReactionExchJoules < 0.0)
        {
            activationEnergy -= heatOfReactionExchJoules;
        }

        
        
        if(EcPQ > activationEnergy) // i.e. exchange reaction can possibly occur.
        {
            scalar summation = 0.0; // declare "summation" term.

            if(activationEnergy < cloud_.constProps(typeIdP).thetaV()[0]*physicoChemical::k.value())
            {
                summation = 1.0; // this refers to the first sentence in Bird's QK paper after Eq.(12).
            }
            else
            {
                iaP = EcPQ / (physicoChemical::k.value()*cloud_.constProps(typeIdP).thetaV()[0]);

                for(label i = 0 ; i <= iaP ; i++)
                {
                    summation += pow((1.0 - ((i*physicoChemical::k.value()*cloud_.constProps(typeIdP).thetaV()[0])/EcPQ)),(1.5 - omegaPQ));
                }
            }

            P_exch = pow((1.0 - (activationEnergy/EcPQ)),(1.5 - omegaPQ))/summation;// now based on modified activation energy.
            
            totalReactionProbability += P_exch;
            reactionProbabilities[3] = P_exch;
        }
       
        if(chargeExchange_)
        {
            //calculate charge exchange probability
            
            label maxElectronicLevelP = cloud_.constProps(typeIdP).numberOfElectronicLevels();
            
            scalar heatOfReactionChargeExchJoules = heatOfReactionChargeExchange_*physicoChemical::k.value();
            
            scalar EcP = translationalEnergy + EEleP;
            
            scalar aDash = aCoeffCharge_*(pow(2.5 - omegaPQ, bCoeffCharge_)*exp(lgamma(2.5 - omegaPQ))/exp(lgamma(2.5 - omegaPQ + bCoeffCharge_)));
            
            TColl = (translationalEnergy/(physicoChemical::k.value()))/(2.5 - omegaPQ);
            
            scalar activationEnergy =  (aDash*pow((TColl/273.0) , bCoeffCharge_) * fabs(heatOfReactionChargeExchJoules));
            
    //         scalar activationEnergy = activationEnergy_ + (aCoeff_*pow((5000/273.0) , bCoeff_) * fabs(heatOfReactionExchJoules));
            
            if(heatOfReactionChargeExchJoules < 0.0)
            {
                activationEnergy -= heatOfReactionChargeExchJoules;
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
                
                EcP = translationalEnergy + EEleP /*+ heatOfReactionExchJoules*/;
                
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
                    
                    label nState = ceil(cloud_.rndGen().sample01<scalar>()*(nPossStates));
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
                    
                    if(prob > cloud_.rndGen().sample01<scalar>())
                    {
                        //Charge exchange can occur
                        totalReactionProbability += prob;
                        reactionProbabilities[4] = prob;
                    }
                }
            }
            
        }
        
        //Decide if a reaction is to occur
        
        if(totalReactionProbability > cloud_.rndGen().sample01<scalar>())
        {
            //A chemical reaction is to occur, choose which one
            
            scalarList normalisedProbabilities(reactionProbabilities.size(), 0.0);
            scalar cumulativeProbability = 0.0;
            
            forAll(normalisedProbabilities, i)
            {
                normalisedProbabilities[i] = reactionProbabilities[i]/totalReactionProbability;
                
                //If current reaction can't occur, don't check for it
                if(normalisedProbabilities[i] > VSMALL)
                {
                    cumulativeProbability += normalisedProbabilities[i];
                    
                    if(cumulativeProbability > cloud_.rndGen().sample01<scalar>())
                    {
                        //Current reaction is to occur
                        
                        if(i == 0)
                        {
                            //Dissociation is to occur
                            dissocReaction = true;
                            break;
                        }
                        if(i == 1)
                        {
                            //Molecule ionisation reaction is to occur
                            ionisationReactionP = true;
                            break;
                        }
                        if(i == 2)
                        {
                            //Atom ionisation reaction is to occur
                            ionisationReactionQ = true;
                            break;
                        }
                        if(i == 3)
                        {
                            //Exchange reaction is to occur
                            exchangeReaction = true;
                            break;
                        }
                        if(i == 4)
                        {
                            //Exchange reaction is to occur
                            chargeExchangeReaction = true;
                            break;
                        }
                    }
                }
            }
        }
        
        //Perform a dissociation reaction
        if(dissocReaction)
        {
            nTotDissociationReactions_++;
            nDissociationReactionsPerTimeStep_++;
            
            if (allowSplitting_)
            {
                relax_ = false; 
                
                scalar heatOfReactionDissJoules = heatOfReactionDiss_*physicoChemical::k.value();
                
                translationalEnergy = translationalEnergy + heatOfReactionDissJoules + EVibP;
                
                translationalEnergy += EEleQ;
                    
                label ELevelQ = cloud_.postCollisionElectronicEnergyLevel
                                (
                                    translationalEnergy,
                                    jMaxQ,
                                    omegaPQ,
                                    EElistQ,
                                    gListQ
                                );
                                
                translationalEnergy -= EElistQ[ELevelQ];
                
                scalar relVelNonDissoMol = sqrt(2.0*translationalEnergy/mR);

                //center of mass velocity of all particles
                vector Ucm = (mP*UP + mQ*UQ)/(mP + mQ);

                // Variable Hard Sphere collision part

                scalar cosTheta = 2.0*cloud_.rndGen().sample01<scalar>() - 1.0;
            
                scalar sinTheta = sqrt(1.0 - cosTheta*cosTheta);
            
                scalar phi = twoPi*cloud_.rndGen().sample01<scalar>();
            
                vector postCollisionRelU =
                    relVelNonDissoMol
                    *vector
                        (
                            cosTheta,
                            sinTheta*cos(phi),
                            sinTheta*sin(phi)
                        );

                //P dissociates
                UP = Ucm + postCollisionRelU*mQ/(mP + mQ);
                UQ = Ucm - postCollisionRelU*mP/(mP + mQ);

                const label& typeId1 = dissociationProductIds_[0];
                const label& typeId2 = dissociationProductIds_[1];

                //Mass of Product one and two
                scalar mP1 = cloud_.constProps(typeId1).mass();
                scalar mP2 = cloud_.constProps(typeId2).mass();

                scalar mRatoms = mP1*mP2/(mP1 + mP2);

                //center of mass velocity of all particles

                vector UcmAtoms = UP;
                
                scalar cRatoms = sqrt(2.0*(ERotP+EEleP)/mRatoms);

                // Variable Hard Sphere collision part
            
                scalar cosTheta2 = 2.0*cloud_.rndGen().sample01<scalar>() - 1.0;
            
                scalar sinTheta2 = sqrt(1.0 - cosTheta2*cosTheta2);
            
                scalar phi2 = twoPi*cloud_.rndGen().sample01<scalar>();
            
                vector postCollisionRelU2 = cRatoms
                *vector
                    (
                        cosTheta2,
                        sinTheta2*cos(phi2),
                        sinTheta2*sin(phi2)
                    );

                vector uP1 = UcmAtoms + postCollisionRelU2*mP2/(mP1 + mP2);
                vector uP2 = UcmAtoms - postCollisionRelU2*mP1/(mP1 + mP2);

                // New atom Q velocity.
                q.U() = UQ;
                q.ELevel() = ELevelQ;

                // Molecule P will dissociate into 2 atoms.
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
                p.ERot() = 0.0;
                p.vibLevel().setSize(0,0);
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
            nTotIonisationReactionsP_++;
            nIonisationReactionsPPerTimeStep_++;
            
            if(allowSplitting_)
            {
                relax_ = false;
                
                scalar heatOfReactionIonisationJoules = heatOfReactionIonP_*physicoChemical::k.value();
                
                translationalEnergy = translationalEnergy + heatOfReactionIonisationJoules + EEleP;
                
                translationalEnergy += EEleQ;
                    
                label ELevelQ = cloud_.postCollisionElectronicEnergyLevel
                                (
                                    translationalEnergy,
                                    jMaxQ,
                                    omegaPQ,
                                    EElistQ,
                                    gListQ
                                );
                                
                translationalEnergy -= EElistQ[ELevelQ];
                
                scalar relVelNonDissoMol = sqrt(2.0*translationalEnergy/mR);

                //center of mass velocity of all particles
                vector Ucm = (mP*UP + mQ*UQ)/(mP + mQ);    

                // Variable Hard Sphere collision part

                scalar cosTheta = 2.0*cloud_.rndGen().sample01<scalar>() - 1.0;
            
                scalar sinTheta = sqrt(1.0 - cosTheta*cosTheta);
            
                scalar phi = twoPi*cloud_.rndGen().sample01<scalar>();
            
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

                const label& typeId1 = ionisationProductsIdsP_[0];
                const label& typeId2 = ionisationProductsIdsP_[1];
                
                //Mass of Product one and two
                scalar mP1 = cloud_.constProps(typeId1).mass();
                scalar mP2 = cloud_.constProps(typeId2).mass();
                
                scalar mRatoms = mP1*mP2/(mP1 + mP2);
                
                translationalEnergy = ERotP + EVibP;
                
                scalar cRatoms = sqrt(2.0*translationalEnergy/mRatoms);

                // Variable Hard Sphere collision part
                scalar cosTheta2 = 2.0*cloud_.rndGen().sample01<scalar>() - 1.0;
            
                scalar sinTheta2 = sqrt(1.0 - cosTheta2*cosTheta2);
            
                scalar phi2 = twoPi*cloud_.rndGen().sample01<scalar>();
            
                vector postCollisionRelU2 = cRatoms
                *vector
                    (
                        cosTheta2,
                        sinTheta2*cos(phi2),
                        sinTheta2*sin(phi2)
                    );


                vector uP1 = UP + postCollisionRelU2*mP2/(mP1 + mP2);
                vector uP2 = UP - postCollisionRelU2*mP1/(mP1 + mP2);

                // Q remains NON-IONISED.
                q.U() = UQ;
                q.ELevel() = ELevelQ;

                // Molecule P will dissociation.
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
                p.vibLevel()[0] = 0;
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
        
        if(ionisationReactionQ)
        {
           nTotIonisationReactionsQ_++;
           nIonisationReactionsQPerTimeStep_++;
           
           if(allowSplitting_)
            {
                relax_ = false;
                
                scalar heatOfReactionIonisationJoules = heatOfReactionIonQ_*physicoChemical::k.value();
                
                translationalEnergy = translationalEnergy + heatOfReactionIonisationJoules + EEleQ;
                
                translationalEnergy += EEleP;
                    
                label ELevelP = cloud_.postCollisionElectronicEnergyLevel
                                (
                                    translationalEnergy,
                                    jMaxP,
                                    omegaPQ,
                                    EElistP,
                                    gListP
                                );
                                
                translationalEnergy -= EElistP[ELevelP];
                
                translationalEnergy += EVibP;
                
                label iMax = (translationalEnergy / (physicoChemical::k.value()*thetaVP));
                
                label vibLevelP = cloud_.postCollisionVibrationalEnergyLevel
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
                
                scalar relVelNonDissoMol = sqrt(2.0*translationalEnergy/mR);

                //center of mass velocity of all particles
                vector Ucm = (mP*UP + mQ*UQ)/(mP + mQ);    

                // Variable Hard Sphere collision part

                scalar cosTheta = 2.0*cloud_.rndGen().sample01<scalar>() - 1.0;
            
                scalar sinTheta = sqrt(1.0 - cosTheta*cosTheta);
            
                scalar phi = twoPi*cloud_.rndGen().sample01<scalar>();
            
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

                const label& typeId1 = ionisationProductsIdsQ_[0];
                const label& typeId2 = ionisationProductsIdsQ_[1];
                
                //Mass of Product one and two
                scalar mP1 = cloud_.constProps(typeId1).mass();
                scalar mP2 = cloud_.constProps(typeId2).mass();
                
                scalar mRatoms = mP1*mP2/(mP1 + mP2);
                
                translationalEnergy = 0; //ERotQ + EVibQ
                
                scalar cRatoms = sqrt(2.0*translationalEnergy/mRatoms);

                // Variable Hard Sphere collision part
                scalar cosTheta2 = 2.0*cloud_.rndGen().sample01<scalar>() - 1.0;
            
                scalar sinTheta2 = sqrt(1.0 - cosTheta2*cosTheta2);
            
                scalar phi2 = twoPi*cloud_.rndGen().sample01<scalar>();
            
                vector postCollisionRelU2 = cRatoms
                *vector
                    (
                        cosTheta2,
                        sinTheta2*cos(phi2),
                        sinTheta2*sin(phi2)
                    );


                vector uQ1 = UQ + postCollisionRelU2*mP2/(mP1 + mP2);
                vector uQ2 = UQ - postCollisionRelU2*mP1/(mP1 + mP2);

                // P remains NON-IONISED.
                p.U() = UP;
                p.ERot() = ERotP;
                p.vibLevel() = vibLevelP;
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
                    typeId2,
                    -1,
                    classificationQ,
                    vibLevel
                );
            }
        }
        
        //Perform exchange reaction
        if(exchangeReaction)
        {
            nTotExchangeReactions_++;
            nExchangeReactionsPerTimeStep_++;
            
            if (allowSplitting_)
            {
                relax_ = false;
                
                //center of mass velocity (pre-exchange)
                vector Ucm = (mP*UP + mQ*UQ)/(mP + mQ);

                const label& typeIdMol = exchangeProductIds_[0];
                const label& typeIdAtom = exchangeProductIds_[1];

                // change species properties

                scalar mPExch = cloud_.constProps(typeIdAtom).mass();
                scalar mQExch = cloud_.constProps(typeIdMol).mass();

                scalar mRExch = mPExch*mQExch/(mPExch + mQExch);
                
                translationalEnergy = translationalEnergy + ERotP + EVibP + EEleP + EEleQ + heatOfReactionExchJoules;
                
                scalar relVelExchMol = sqrt((2.0*translationalEnergy)/mRExch);
                
                //Variable Hard Sphere collision part for collision of molecules
        
                scalar cosTheta = 2.0*cloud_.rndGen().sample01<scalar>() - 1.0;
            
                scalar sinTheta = sqrt(1.0 - cosTheta*cosTheta);
            
                scalar phi = twoPi*cloud_.rndGen().sample01<scalar>();
            
                vector postCollisionRelU =
                    relVelExchMol
                    *vector
                        (
                            cosTheta,
                            sinTheta*cos(phi),
                            sinTheta*sin(phi)
                        );
                
                q.typeId() = typeIdMol; // q is originally the atom, becomes the molecule
                p.typeId() = typeIdAtom; // p is the originally the molecule, becomes the atom
        
                UP = Ucm + (postCollisionRelU*mQExch/(mPExch + mQExch)); // P changes from mol to atom.
                UQ = Ucm - (postCollisionRelU*mPExch/(mPExch + mQExch)); // Q changes from atom to mol.

                if(
                    cloud_.constProps(exchangeProductIds_[0]).
                    rotationalDegreesOfFreedom() > VSMALL
                )
                {
                    q.vibLevel().setSize(1,0);
                }       
                else
                {
                    q.vibLevel().setSize(0,0);
                }
               
                q.ELevel() = 0;
                q.ERot() = 0;
                q.U() = UQ;

                p.U() = UP;
                p.ERot() = 0.0; // remove p's internal energies as it's now an atom
                if(
                    cloud_.constProps(exchangeProductIds_[1]).
                    rotationalDegreesOfFreedom() > VSMALL
                )
                {
                    p.vibLevel().setSize(1,0);
                }       
                else
                {
                    p.vibLevel().setSize(0,0);
                }
                p.ELevel() = 0;
            }   
        }
        
        if(chargeExchangeReaction)
        {
            relax_ = false;
            
            nTotChargeExchangeReactions_++;
            nChargeExchangeReactionsPerTimeStep_++;
            
            if(allowSplitting_)
            {
                relax_ = false;
                
                scalar heatOfReactionChargeExchJoules = heatOfReactionChargeExchange_*physicoChemical::k.value();
                
                translationalEnergy = translationalEnergy + heatOfReactionChargeExchJoules;
                
                translationalEnergy += ERotP + EVibP + EEleP + ERotQ  + EEleQ;
                //+EVibQ
                scalar relVel = sqrt((2.0*translationalEnergy)/mR);
                    
                // centre of mass velocity of molecules (pre-split)
                vector Ucm = (mP*UP + mQ*UQ)/(mP + mQ);

                // Variable Hard Sphere collision part for collision of molecules
                scalar cosTheta = 2.0*cloud_.rndGen().sample01<scalar>() - 1.0;
            
                scalar sinTheta = sqrt(1.0 - cosTheta*cosTheta);
            
                scalar phi = twoPi*cloud_.rndGen().sample01<scalar>();
            
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
             
    // if q is the molecule and p is the atom... 
    if(typeIdP == reactantIds_[1] && typeIdQ == reactantIds_[0])// C + AB ---> A + B + C (diss type II) or C + AB ---> AC + B (forward exchange)
    {
        relax_ = true;
        
        scalar totalReactionProbability = 0.0;
        scalarList reactionProbabilities(5, 0.0);
        
        vector UP = p.U();
        vector UQ = q.U();
        
        scalar ERotP = p.ERot();
        scalar ERotQ = q.ERot();
        
        
        scalar EEleP = cloud_.constProps(typeIdP).electronicEnergyList()[p.ELevel()];
        scalar EEleQ = 
            cloud_.constProps(typeIdQ).electronicEnergyList()[q.ELevel()];
       
        scalar EVibQ = 
                q.vibLevel()[0]*(cloud_.constProps(typeIdQ).thetaV()[0]
                *physicoChemical::k.value());
        
        scalar mP = cloud_.constProps(typeIdP).mass();
        scalar mQ = cloud_.constProps(typeIdQ).mass();
        
        label jMaxP = cloud_.constProps(typeIdP).numberOfElectronicLevels();
        List<label> gListP = cloud_.constProps(typeIdP).degeneracyList();
        List<scalar> EElistP = cloud_.constProps(typeIdP).electronicEnergyList();
        
        label jMaxQ = cloud_.constProps(typeIdQ).numberOfElectronicLevels();
        List<label> gListQ = cloud_.constProps(typeIdQ).degeneracyList();
        List<scalar> EElistQ = cloud_.constProps(typeIdQ).electronicEnergyList();

        scalar mR = mP*mQ/(mP + mQ);
        scalar cRsqr = magSqr(UP - UQ);
        scalar translationalEnergy = 0.5*mR*cRsqr;
        
        scalar thetaVQ = cloud_.constProps(typeIdQ).thetaV()[0];
        scalar thetaDQ = cloud_.constProps(typeIdQ).thetaD()[0];
        scalar ZrefQ = cloud_.constProps(typeIdQ).Zref()[0];
        scalar refTempZvQ = cloud_.constProps(typeIdQ).TrefZv()[0];
        scalar rotationalDofQ = cloud_.constProps(typeIdQ).rotationalDegreesOfFreedom();
        
        scalar omegaPQ =
        0.5
        *(
            cloud_.constProps(typeIdP).omega()
            + cloud_.constProps(typeIdQ).omega()
        );
        
        scalar ChiB = 2.5 - omegaPQ;

        scalar EcPQ = 0.0;
        scalar TColl = 0.0;
        label idQ = -1;
        label deltaDissoIQ= 0;
        label imaxQ = 0;
        label iaQ = 0;
        bool dissocReaction = false;
        bool ionisationReactionP = false;
        bool ionisationReactionQ = false;
        bool exchangeReaction = false;
        bool chargeExchangeReaction = false;
        
        
        if(!chargedMolecule_)
        {
            // firstly calculate dissociation probability (0 or 1).
            EcPQ = translationalEnergy + EVibQ;

            imaxQ = EcPQ/(physicoChemical::k.value()*cloud_.constProps(typeIdQ).thetaV()[0]);
            
            idQ = cloud_.constProps(typeIdQ).thetaD()[0]/cloud_.constProps(typeIdQ).thetaV()[0];
            
            deltaDissoIQ = imaxQ-idQ;
            
            if (deltaDissoIQ > 0)
            {
                // DISSOCIATION CAN OCCUR//
                totalReactionProbability += 1.0;
                reactionProbabilities[0] = 1.0;
            }
            
            //Now, ionisation of the molecule (Q)
            
            scalar ionisationEnergy = cloud_.constProps(typeIdQ).ionisationTemperature()*physicoChemical::k.value();
            
            // calculate if an ionisation of species Q is possible
            EcPQ = translationalEnergy + EEleQ;

            if((EcPQ - ionisationEnergy) > VSMALL)
            {
                //Ionisation can occur
                totalReactionProbability += 1.0;
                reactionProbabilities[1] = 1.0;
            }
        }
        
        if(!chargedAtom_)
        {
            //Now, ionisation of the atom (P)
            
            scalar ionisationEnergy = cloud_.constProps(typeIdP).ionisationTemperature()*physicoChemical::k.value();
            
            // calculate if an ionisation of species Q is possible
            EcPQ = translationalEnergy + EEleP;

            if((EcPQ - ionisationEnergy) > VSMALL)
            {
                //Ionisation can occur
                totalReactionProbability += 1.0;
                reactionProbabilities[2] = 1.0;
            }
        }     
        
        EcPQ = translationalEnergy + EVibQ;
        
        scalar P_exch = 0.0;
            
        TColl = (translationalEnergy/(physicoChemical::k.value()))/(2.5 - omegaPQ);// Eq.(10) of Bird QK paper.
            
        scalar heatOfReactionExchJoules = heatOfReactionExch_*physicoChemical::k.value();
        
        scalar aDash = aCoeff_*(pow(2.5 - omegaPQ, bCoeff_)*exp(lgamma(2.5 - omegaPQ))/exp(lgamma(2.5 - omegaPQ + bCoeff_)));

        scalar activationEnergy = (aDash*pow((TColl/273.0) , bCoeff_) * fabs(heatOfReactionExchJoules));
    
        if(heatOfReactionExchJoules < 0.0)
        {
            activationEnergy -= heatOfReactionExchJoules;
        }
        
        if(EcPQ > activationEnergy) // i.e. exchange reaction can possibly occur.
        {
            scalar summation = 0.0; // declare "summation" term.

            if(activationEnergy < cloud_.constProps(typeIdQ).thetaV()[0]*physicoChemical::k.value())
            {
                summation = 1.0; // this refers to the first sentence in Bird's QK paper after Eq.(12).
            }
            else
            {
                iaQ = EcPQ / (physicoChemical::k.value()*cloud_.constProps(typeIdQ).thetaV()[0]);

                for(label i = 0 ; i <= iaQ ; i++)
                {
                    summation += pow((1.0 - ((i*physicoChemical::k.value()*cloud_.constProps(typeIdQ).thetaV()[0])/EcPQ)),(1.5 - omegaPQ));
                }
            }
        
            P_exch = pow((1.0 - (activationEnergy/EcPQ)),(1.5 - omegaPQ))/summation;// now based on modified activation energy.
            
            totalReactionProbability += P_exch;
            reactionProbabilities[3] = P_exch;
        }
        
        
        
        if(chargeExchange_)
        {
            label maxElectronicLevelQ = cloud_.constProps(typeIdQ).numberOfElectronicLevels();
            
            scalar heatOfReactionChargeExchJoules = heatOfReactionChargeExchange_*physicoChemical::k.value();
        
            scalar EcQ = translationalEnergy + EEleQ;
            
            scalar aDash = aCoeffCharge_*(pow(2.5 - omegaPQ, bCoeffCharge_)*exp(lgamma(2.5 - omegaPQ))/exp(lgamma(2.5 - omegaPQ + bCoeffCharge_)));
            
            TColl = (translationalEnergy/(physicoChemical::k.value()))/(2.5 - omegaPQ);

            scalar activationEnergy = (aDash*pow((TColl/273.0) , bCoeffCharge_) * fabs(heatOfReactionChargeExchJoules));
            
    //         scalar activationEnergy = (aCoeff_*pow((5000/273.0) , bCoeff_) * fabs(heatOfReactionExchJoules));

            if(heatOfReactionChargeExchJoules < 0.0)
            {
                activationEnergy -= heatOfReactionChargeExchJoules;
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
                
                EcQ = translationalEnergy + EEleQ /*+ heatOfReactionExchJoules*/;

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
                    
                    label nState = ceil(cloud_.rndGen().sample01<scalar>()*(nPossStates));
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
                    
                    if(prob > cloud_.rndGen().sample01<scalar>())
                    {
                        //Charge exchange can occur
                        totalReactionProbability += prob;
                        reactionProbabilities[4] = prob;
                    }
                }
            }
        }
               
        //Decide if a reaction is to occur
        
        if(totalReactionProbability > cloud_.rndGen().sample01<scalar>())
        {
            //A chemical reaction is to occur, choose which one
            
            scalarList normalisedProbabilities(reactionProbabilities.size(), 0.0);
            scalar cumulativeProbability = 0.0;
            
            forAll(normalisedProbabilities, i)
            {
                normalisedProbabilities[i] = reactionProbabilities[i]/totalReactionProbability;
                
                //If current reaction can't occur, don't check for it
                if(normalisedProbabilities[i] > VSMALL)
                {
                    cumulativeProbability += normalisedProbabilities[i];
                    
                    if(cumulativeProbability > cloud_.rndGen().sample01<scalar>())
                    {
                        //Current reaction is to occur
                        
                        if(i == 0)
                        {
                            //Dissociation is to occur
                            dissocReaction = true;
                            break;
                        }
                        if(i == 1)
                        {
                            //Molecule ionisation reaction is to occur
                            ionisationReactionP = true;
                            break;
                        }
                        if(i == 2)
                        {
                            //Atom ionisation reaction is to occur
                            ionisationReactionQ = true;
                            break;
                        }
                        if(i == 3)
                        {
                            //Exchange reaction is to occur
                            exchangeReaction = true;
                            break;
                        }
                        if(i == 4)
                        {
                            //Exchange reaction is to occur
                            chargeExchangeReaction = true;
                            break;
                        }
                    }
                }
            }
        }
        
        //Perform dissociation reaction
        if(dissocReaction)
        {
            nTotDissociationReactions_++;
            nDissociationReactionsPerTimeStep_++;
            
            if (allowSplitting_)
            {
                relax_ = false; 
                
                scalar heatOfReactionDissJoules = heatOfReactionDiss_*physicoChemical::k.value();
                
                translationalEnergy = translationalEnergy + heatOfReactionDissJoules + EVibQ;
                
                translationalEnergy += EEleP;
                    
                label ELevelP = cloud_.postCollisionElectronicEnergyLevel
                                (
                                    translationalEnergy,
                                    jMaxP,
                                    omegaPQ,
                                    EElistP,
                                    gListP
                                );
                                
                translationalEnergy -= EElistP[ELevelP];
                
                scalar relVelNonDissoMol = sqrt(2.0*translationalEnergy/mR);

                //center of mass velocity of all particles
                vector Ucm = (mP*UP + mQ*UQ)/(mP + mQ);

                // Variable Hard Sphere collision part

                scalar cosTheta = 2.0*cloud_.rndGen().sample01<scalar>() - 1.0;
            
                scalar sinTheta = sqrt(1.0 - cosTheta*cosTheta);
            
                scalar phi = twoPi*cloud_.rndGen().sample01<scalar>();
            
                vector postCollisionRelU =
                    relVelNonDissoMol
                    *vector
                        (
                            cosTheta,
                            sinTheta*cos(phi),
                            sinTheta*sin(phi)
                        );


                UP = Ucm + postCollisionRelU*mQ/(mP + mQ); // P is the NON-DISSOCIATING molecule.
                UQ = Ucm - postCollisionRelU*mP/(mP + mQ); 

                const label& typeId1 = dissociationProductIds_[0];
                const label& typeId2 = dissociationProductIds_[1];

                //Mass of Product one and two
                scalar mP1 = cloud_.constProps(typeId1).mass();
                scalar mP2 = cloud_.constProps(typeId2).mass();

                scalar mRatoms = mP1*mP2/(mP1 + mP2);

                //center of mass velocity of all particles

                vector UcmAtoms = UQ;
                
                scalar cRatoms = sqrt(2.0*(ERotQ+EEleQ)/mRatoms);

                // Variable Hard Sphere collision part
            
                scalar cosTheta2 = 2.0*cloud_.rndGen().sample01<scalar>() - 1.0;
            
                scalar sinTheta2 = sqrt(1.0 - cosTheta2*cosTheta2);
            
                scalar phi2 = twoPi*cloud_.rndGen().sample01<scalar>();
            
                vector postCollisionRelU2 = cRatoms
                *vector
                    (
                        cosTheta2,
                        sinTheta2*cos(phi2),
                        sinTheta2*sin(phi2)
                    );

                vector uQ1 = UcmAtoms + postCollisionRelU2*mP2/(mP1 + mP2);
                vector uQ2 = UcmAtoms - postCollisionRelU2*mP1/(mP1 + mP2);

                // New atom P velocity.
                p.U() = UP;
                p.ELevel() = ELevelP;

                // Molecule Q will dissociate into 2 atoms.
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
                q.ERot() = 0.0;
                q.vibLevel().setSize(0,0);
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
        
        if(ionisationReactionP)
        {
            //Molecule ionisation (Q is the molecule, P is used for measurement purposes)
            nTotIonisationReactionsP_++;
            nIonisationReactionsPPerTimeStep_++;
            
            if(allowSplitting_)
            {
                relax_ = false;
                
                scalar heatOfReactionIonisationJoules = heatOfReactionIonP_*physicoChemical::k.value();
                
                translationalEnergy = translationalEnergy + heatOfReactionIonisationJoules + EEleQ;
                
                translationalEnergy += EEleP;
                    
                label ELevelP = cloud_.postCollisionElectronicEnergyLevel
                                (
                                    translationalEnergy,
                                    jMaxP,
                                    omegaPQ,
                                    EElistP,
                                    gListP
                                );
                                
                translationalEnergy -= EElistP[ELevelP];
                
                scalar relVelNonDissoMol = sqrt(2.0*translationalEnergy/mR);

                //center of mass velocity of all particles
                vector Ucm = (mP*UP + mQ*UQ)/(mP + mQ);    

                // Variable Hard Sphere collision part

                scalar cosTheta = 2.0*cloud_.rndGen().sample01<scalar>() - 1.0;
            
                scalar sinTheta = sqrt(1.0 - cosTheta*cosTheta);
            
                scalar phi = twoPi*cloud_.rndGen().sample01<scalar>();
            
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

                const label& typeId1 = ionisationProductsIdsP_[0];
                const label& typeId2 = ionisationProductsIdsP_[1];
                
                //Mass of Product one and two
                scalar mP1 = cloud_.constProps(typeId1).mass();
                scalar mP2 = cloud_.constProps(typeId2).mass();
                
                scalar mRatoms = mP1*mP2/(mP1 + mP2);
                
                translationalEnergy = ERotQ + EVibQ;
                
                scalar cRatoms = sqrt(2.0*translationalEnergy/mRatoms);

                // Variable Hard Sphere collision part
                scalar cosTheta2 = 2.0*cloud_.rndGen().sample01<scalar>() - 1.0;
            
                scalar sinTheta2 = sqrt(1.0 - cosTheta2*cosTheta2);
            
                scalar phi2 = twoPi*cloud_.rndGen().sample01<scalar>();
            
                vector postCollisionRelU2 = cRatoms
                *vector
                    (
                        cosTheta2,
                        sinTheta2*cos(phi2),
                        sinTheta2*sin(phi2)
                    );


                vector uQ1 = UQ + postCollisionRelU2*mP2/(mP1 + mP2);
                vector uQ2 = UQ - postCollisionRelU2*mP1/(mP1 + mP2);

                // P remains NON-IONISED.
                p.U() = UP;
                p.ELevel() = ELevelP;

                // Molecule Q will ionise
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
                q.vibLevel().setSize(1,0);
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
        
        if(ionisationReactionQ)
        {
           //P is the atom (Q is used for measurement purposes) 
           nTotIonisationReactionsQ_++;
           nIonisationReactionsQPerTimeStep_++;
           
           if(allowSplitting_)
            {
                relax_ = false;
                
                scalar heatOfReactionIonisationJoules = heatOfReactionIonQ_*physicoChemical::k.value();
                
                translationalEnergy = translationalEnergy + heatOfReactionIonisationJoules + EEleP;
                
                translationalEnergy += EEleQ;
                    
                label ELevelQ = cloud_.postCollisionElectronicEnergyLevel
                                (
                                    translationalEnergy,
                                    jMaxQ,
                                    omegaPQ,
                                    EElistQ,
                                    gListQ
                                );
                                
                translationalEnergy -= EElistQ[ELevelQ];
                
                translationalEnergy += EVibQ;
                
                label iMax = (translationalEnergy / (physicoChemical::k.value()*thetaVQ));
                
                label vibLevelQ = cloud_.postCollisionVibrationalEnergyLevel
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
                
                scalar relVelNonDissoMol = sqrt(2.0*translationalEnergy/mR);

                //center of mass velocity of all particles
                vector Ucm = (mP*UP + mQ*UQ)/(mP + mQ);    

                // Variable Hard Sphere collision part

                scalar cosTheta = 2.0*cloud_.rndGen().sample01<scalar>() - 1.0;
            
                scalar sinTheta = sqrt(1.0 - cosTheta*cosTheta);
            
                scalar phi = twoPi*cloud_.rndGen().sample01<scalar>();
            
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

                const label& typeId1 = ionisationProductsIdsQ_[0];
                const label& typeId2 = ionisationProductsIdsQ_[1];
                
                //Mass of Product one and two
                scalar mP1 = cloud_.constProps(typeId1).mass();
                scalar mP2 = cloud_.constProps(typeId2).mass();
                
                scalar mRatoms = mP1*mP2/(mP1 + mP2);
                
                translationalEnergy = 0; //ERotP + EVibP
                
                scalar cRatoms = sqrt(2.0*translationalEnergy/mRatoms);

                // Variable Hard Sphere collision part
                scalar cosTheta2 = 2.0*cloud_.rndGen().sample01<scalar>() - 1.0;
            
                scalar sinTheta2 = sqrt(1.0 - cosTheta2*cosTheta2);
            
                scalar phi2 = twoPi*cloud_.rndGen().sample01<scalar>();
            
                vector postCollisionRelU2 = cRatoms
                *vector
                    (
                        cosTheta2,
                        sinTheta2*cos(phi2),
                        sinTheta2*sin(phi2)
                    );


                vector uP1 = UP + postCollisionRelU2*mP2/(mP1 + mP2);
                vector uP2 = UP - postCollisionRelU2*mP1/(mP1 + mP2);

                // Q remains NON-IONISED.
                q.U() = UQ;
                q.ERot() = ERotQ;
                q.vibLevel()[0] = vibLevelQ;
                q.ELevel() = ELevelQ;

                // Atom P will ionise.
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
        
        //Perform exchange reaction
        if(exchangeReaction)
        {          
            nTotExchangeReactions_++;
            nExchangeReactionsPerTimeStep_++;
            
            if(allowSplitting_)
            {   
                relax_ = false;
                
                //center of mass velocity (pre-exchange)
                vector Ucm = (mP*UP + mQ*UQ)/(mP + mQ);

                const label& typeIdMol = exchangeProductIds_[0];
                const label& typeIdAtom = exchangeProductIds_[1];

                // change species properties

                scalar mQExch = cloud_.constProps(typeIdAtom).mass();
                scalar mPExch = cloud_.constProps(typeIdMol).mass();

                scalar mRExch = mPExch*mQExch/(mPExch + mQExch);
                
                translationalEnergy = translationalEnergy + ERotQ + EVibQ + EEleQ + EEleP + heatOfReactionExchJoules;
                
                p.typeId() = typeIdMol; // p is originally the atom, becomes the molecule
                q.typeId() = typeIdAtom; // q is orinally the molecule, becomes the atom

                scalar relVelExchMol = sqrt((2.0*translationalEnergy)/mRExch);

                // Variable Hard Sphere collision part for collision of molecules
        
                scalar cosTheta = 2.0*cloud_.rndGen().sample01<scalar>() - 1.0;
            
                scalar sinTheta = sqrt(1.0 - cosTheta*cosTheta);
            
                scalar phi = twoPi*cloud_.rndGen().sample01<scalar>();
            
                vector postCollisionRelU =
                    relVelExchMol
                    *vector
                        (
                            cosTheta,
                            sinTheta*cos(phi),
                            sinTheta*sin(phi)
                        );
        
                UP = Ucm + (postCollisionRelU*mQExch/(mPExch + mQExch)); // P changes from atom to mol.
                UQ = Ucm - (postCollisionRelU*mPExch/(mPExch + mQExch)); // Q changes from mol to atom. 
                
                
                if(
                    cloud_.constProps(exchangeProductIds_[0]).
                    rotationalDegreesOfFreedom() > VSMALL
                )
                {
                    p.vibLevel().setSize(1,0);
                }       
                else
                {
                    p.vibLevel().setSize(0,0);
                } 
                p.ELevel() = 0;
                p.ERot() = 0;
                p.U() = UP;
                
                q.U() = UQ;
                q.ELevel() = 0;
                if(
                    cloud_.constProps(exchangeProductIds_[1]).
                    rotationalDegreesOfFreedom() > VSMALL
                )
                {
                    q.vibLevel().setSize(1,0);
                }       
                else
                {
                    q.vibLevel().setSize(0,0);
                } 
            }
        }
        
        if(chargeExchangeReaction)
        {
            nTotChargeExchangeReactions_++;
            nChargeExchangeReactionsPerTimeStep_++;
                    
            if(allowSplitting_)
            {
                relax_ = false;
                
                scalar heatOfReactionChargeExchJoules = heatOfReactionChargeExchange_*physicoChemical::k.value();
                
                translationalEnergy = translationalEnergy + heatOfReactionChargeExchJoules;
                
                translationalEnergy += ERotP + EEleP + ERotQ + EVibQ + EEleQ;
                //+EVibP
                scalar relVel = sqrt((2.0*translationalEnergy)/mR);
                    
                    // centre of mass velocity of molecules (pre-split)
                vector Ucm = (mP*UP + mQ*UQ)/(mP + mQ);

                // Variable Hard Sphere collision part for collision of molecules
                scalar cosTheta = 2.0*cloud_.rndGen().sample01<scalar>() - 1.0;
            
                scalar sinTheta = sqrt(1.0 - cosTheta*cosTheta);
            
                scalar phi = twoPi*cloud_.rndGen().sample01<scalar>();
            
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

void dissociationIonisationExchange::reactExchangeMolecule
(
    dsmcParcel& p,
    label newTypeId,
    const label& newEVibLevel,
    const scalar& newERot,
    const vector& newU
)
{
    p.vibLevel() = newEVibLevel;
    p.ERot() = newERot;
    p.U() = newU;
    p.typeId() = newTypeId;
}

void dissociationIonisationExchange::reactExchangeAtom
(
    dsmcParcel& p,
    label newTypeId,
    const vector& newU
)
{
    p.U() = newU;
    p.typeId() = newTypeId;
}

void  dissociationIonisationExchange::outputResults(const label& counterIndex)
{
    if(writeRatesToTerminal_ == true)
    {
        
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
        label nTotExchangeReactions = nTotExchangeReactions_;
        label nTotChargeExchangeReactions = nTotChargeExchangeReactions_;
        label nTotDissociationReactions = nTotDissociationReactions_;
        label nTotIonisationReactionsP = nTotIonisationReactionsP_;
        label nTotIonisationReactionsQ = nTotIonisationReactionsQ_; 

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
            reduce(nTotExchangeReactions, sumOp<label>());
            reduce(nTotChargeExchangeReactions, sumOp<label>());
            reduce(nTotDissociationReactions, sumOp<label>());
            reduce(nTotIonisationReactionsP, sumOp<label>());
            reduce(nTotIonisationReactionsQ, sumOp<label>());
        }

        numberDensities_[0] = (mols[0]*cloud().nParticle())/volume;
        numberDensities_[1] = (mols[1]*cloud().nParticle())/volume;

        const scalar& deltaT = mesh_.time().deltaT().value();

        word reactantMolA = cloud_.typeIdList()[reactantIds_[0]];
        word reactantMolB = cloud_.typeIdList()[reactantIds_[1]];

        word exchangeProductMolA = cloud_.typeIdList()[exchangeProductIds_[0]];
        word exchangeProductMolB = cloud_.typeIdList()[exchangeProductIds_[1]];
        
        word dissociationProductMolA;
        word dissociationProductMolB;
        
        word moleculeIonisationProductMolA;
        word moleculeIonisationProductMolB;
        
        word atomIonisationProductMolA;
        word atomIonisationProductMolB;
        
        word chargeExchangeProductMolA;
        word chargeExchangeProductMolB;

        if(!chargedMolecule_)
        {
            dissociationProductMolA = cloud_.typeIdList()[dissociationProductIds_[0]];
            dissociationProductMolB = cloud_.typeIdList()[dissociationProductIds_[1]];
            
            moleculeIonisationProductMolA = cloud_.typeIdList()[ionisationProductsIdsP_[0]];
            moleculeIonisationProductMolB = cloud_.typeIdList()[ionisationProductsIdsP_[1]];
        }
        
        if(!chargedAtom_)
        {
            atomIonisationProductMolA = cloud_.typeIdList()[ionisationProductsIdsQ_[0]];
            atomIonisationProductMolB = cloud_.typeIdList()[ionisationProductsIdsQ_[1]];
        }
        
        if(chargeExchange_)
        {
            chargeExchangeProductMolA = cloud_.typeIdList()[chargeExchangeProductIds_[0]];
            chargeExchangeProductMolB = cloud_.typeIdList()[chargeExchangeProductIds_[1]];
        }
            
        if((numberDensities_[0] > 0.0) && (numberDensities_[1] > 0.0))
        {   

            scalar reactionRate1 = 0.0;
            scalar reactionRate2 = 0.0;
            scalar reactionRate3 = 0.0;
            scalar reactionRate4 = 0.0;
            scalar reactionRate5 = 0.0;

            reactionRate1 =
            (
                nTotExchangeReactions
                * cloud_.nParticle()
            )/(counterIndex*deltaT*numberDensities_[0]* numberDensities_[1]*volume);

            Info << "Exchange reaction " 
                    <<  reactantMolA << " + " << reactantMolB << " --> " 
                    << exchangeProductMolA << " + " << exchangeProductMolB 
                    << ", reaction rate = " << reactionRate1 
                    << endl;
            
            
            if(!chargedMolecule_)
            {
                reactionRate2 =
                (
                    nTotDissociationReactions
                    * cloud_.nParticle()
                    )/(counterIndex*deltaT*numberDensities_[0]* numberDensities_[1]*volume);
                   
                Info << "Type II dissociation reaction " 
                    <<  reactantMolA << " + " << reactantMolB << " --> " 
                    << dissociationProductMolA << " + " << dissociationProductMolB 
                    <<  " + " << reactantMolB 
                    << ", reaction rate = " << reactionRate2
                    << endl;    
                    
                    
                reactionRate3 =
                (
                    nTotIonisationReactionsP
                    * cloud_.nParticle()
                    )/(counterIndex*deltaT*numberDensities_[0]* numberDensities_[1]*volume);
            
                Info << "Ionisation reaction " 
                    <<  reactantMolA << " + " << reactantMolB << " --> " 
                    << moleculeIonisationProductMolA << " + " << moleculeIonisationProductMolB 
                    <<  " + " << reactantMolB 
                    << ", reaction rate = " << reactionRate3
                    << endl;
            }
            
            if(!chargedAtom_)
            {
                reactionRate4 =
                (
                    nTotIonisationReactionsQ
                    * cloud_.nParticle()
                    )/(counterIndex*deltaT*numberDensities_[0]* numberDensities_[1]*volume);
            
                Info << "Ionisation reaction " 
                    <<  reactantMolA << " + " << reactantMolB << " --> " 
                    <<  reactantMolA << " + "
                    << atomIonisationProductMolA << " + " << atomIonisationProductMolB
                    << ", reaction rate = " << reactionRate4
                    << endl;
                    
            }
            
            if(chargeExchange_)
            {
                reactionRate5 =
                (
                    nTotChargeExchangeReactions
                    * cloud_.nParticle()
                    )/(counterIndex*deltaT*numberDensities_[0]* numberDensities_[1]*volume);
            
               Info << "Charge exchange reaction " 
                    <<  reactantMolA << " + " << reactantMolB << " --> "
                    << chargeExchangeProductMolA << " + " << chargeExchangeProductMolB
                    << ", reaction rate = " << reactionRate5
                    << endl;
            }
        }
    }
    else
    {
        word reactantMolA = cloud_.typeIdList()[reactantIds_[0]];
        word reactantMolB = cloud_.typeIdList()[reactantIds_[1]];

        word exchangeProductMolA = cloud_.typeIdList()[exchangeProductIds_[0]];
        word exchangeProductMolB = cloud_.typeIdList()[exchangeProductIds_[1]];
        
        word dissociationProductMolA;
        word dissociationProductMolB;
            
        word moleculeIonisationProductMolA;
        word moleculeIonisationProductMolB;

        if(!chargedMolecule_)
        {
            dissociationProductMolA = cloud_.typeIdList()[dissociationProductIds_[0]];
            dissociationProductMolB = cloud_.typeIdList()[dissociationProductIds_[1]];
            
            moleculeIonisationProductMolA = cloud_.typeIdList()[ionisationProductsIdsP_[0]];
            moleculeIonisationProductMolB = cloud_.typeIdList()[ionisationProductsIdsP_[1]];
        }
        
        word atomIonisationProductMolA;
        word atomIonisationProductMolB;
        
        if(!chargedAtom_)
        {
            atomIonisationProductMolA = cloud_.typeIdList()[ionisationProductsIdsQ_[0]];
            atomIonisationProductMolB = cloud_.typeIdList()[ionisationProductsIdsQ_[1]];
        }
        
        word chargeExchangeProductMolA;
        word chargeExchangeProductMolB;
        
        if(chargeExchange_)
        {
            chargeExchangeProductMolA = cloud_.typeIdList()[chargeExchangeProductIds_[0]];
            chargeExchangeProductMolB = cloud_.typeIdList()[chargeExchangeProductIds_[1]];
        }
            
        label nTotExchangeReactions = nTotExchangeReactions_;
        label nTotChargeExchangeReactions = nTotChargeExchangeReactions_;
        label nTotDissociationReactions = nTotDissociationReactions_;
        label nTotIonisationReactionsP = nTotIonisationReactionsP_;
        label nTotIonisationReactionsQ = nTotIonisationReactionsQ_;
        
        label nDissociationReactionsPerTimeStep = nDissociationReactionsPerTimeStep_;
        label nIonisationReactionsPPerTimeStep = nIonisationReactionsPPerTimeStep_;
        label nIonisationReactionsQPerTimeStep = nIonisationReactionsQPerTimeStep_;
        label nExchangeReactionsPerTimeStep = nExchangeReactionsPerTimeStep_;
        label nChargeExchangeReactionsPerTimeStep = nChargeExchangeReactionsPerTimeStep_;
        
        //- Parallel communication
        if(Pstream::parRun())
        {
            reduce(nTotExchangeReactions, sumOp<label>());
            reduce(nTotChargeExchangeReactions, sumOp<label>());
            reduce(nTotDissociationReactions, sumOp<label>());
            reduce(nTotIonisationReactionsP, sumOp<label>());
            reduce(nTotIonisationReactionsQ, sumOp<label>());
            
            reduce(nDissociationReactionsPerTimeStep, sumOp<label>());
            reduce(nIonisationReactionsPPerTimeStep, sumOp<label>());
            reduce(nIonisationReactionsQPerTimeStep, sumOp<label>());
            reduce(nExchangeReactionsPerTimeStep, sumOp<label>());
            reduce(nChargeExchangeReactionsPerTimeStep, sumOp<label>());
        }
        
        if(nTotExchangeReactions > VSMALL)
        {
            Info << "Exchange reaction " 
                    <<  reactantMolA << " + " << reactantMolB << " --> " 
                    << exchangeProductMolA << " + " << exchangeProductMolB 
                    << " is active, nReactions this time step = " << nExchangeReactionsPerTimeStep << endl;
        }
        
        if(nTotDissociationReactions > VSMALL)
        {
            Info  << "Type II dissociation reaction " 
                <<  reactantMolA << " + " << reactantMolB << " --> " 
                << dissociationProductMolA << " + " << dissociationProductMolB 
                <<  " + " << reactantMolB 
                << " is active, nReactions this time step = " << nDissociationReactionsPerTimeStep << endl;
        } 
        
        if(nTotIonisationReactionsP > VSMALL)
        {
            Info  << "Ionisation reaction " 
                <<  reactantMolA << " + " << reactantMolB << " --> " 
                << moleculeIonisationProductMolA << " + " << moleculeIonisationProductMolB 
                <<  " + " << reactantMolB 
                << " is active, nReactions this time step = " << nIonisationReactionsPPerTimeStep << endl;
        }
        
        if(nTotIonisationReactionsQ > VSMALL)
        {
            Info  << "Ionisation reaction " 
                <<  reactantMolA << " + " << reactantMolB << " --> " 
                <<  reactantMolA << " + "
                << atomIonisationProductMolA << " + " << atomIonisationProductMolB 
                << " is active, nReactions this time step = " << nIonisationReactionsQPerTimeStep << endl;
        }
        
        if(nTotChargeExchangeReactions > VSMALL)
        {
            Info  << "Charge exchange reaction " 
                <<  reactantMolA << " + " << reactantMolB << " --> " 
                << chargeExchangeProductMolA << " + " << chargeExchangeProductMolB 
                << " is active, nReactions this time step = " << nChargeExchangeReactionsPerTimeStep << endl;
        }
    }

    nDissociationReactionsPerTimeStep_ = 0.0;
    nIonisationReactionsPPerTimeStep_ = 0.0;
    nIonisationReactionsQPerTimeStep_ = 0.0;
    nExchangeReactionsPerTimeStep_ = 0.0;
    nChargeExchangeReactionsPerTimeStep_ = 0.0;
}


const bool& dissociationIonisationExchange::relax() const
{
    return relax_;
}

}
// End namespace Foam

// ************************************************************************* //
