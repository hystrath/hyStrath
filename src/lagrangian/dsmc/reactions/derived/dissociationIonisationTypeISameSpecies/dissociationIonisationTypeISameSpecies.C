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

#include "dissociationIonisationTypeISameSpecies.H"
#include "addToRunTimeSelectionTable.H"
#include "fvc.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(dissociationIonisationTypeISameSpecies, 0);

addToRunTimeSelectionTable(dsmcReaction, dissociationIonisationTypeISameSpecies, dictionary);



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
dissociationIonisationTypeISameSpecies::dissociationIonisationTypeISameSpecies
(
    Time& t,
    dsmcCloud& cloud,
    const dictionary& dict
)
:
    dsmcReaction(t, cloud, dict),
    propsDict_(dict.subDict(typeName + "Properties")),
    reactantIds_(),
    productIdsIonisation_(),
    productIdsDissociation_(),
    reactionName_(propsDict_.lookup("reactionName")),
    nTotIonisationReactions_(0),
    nIonisationReactionsPerTimeStep_(0),
    nTotDissociationReactions_(0),
    nDissociationReactionsPerTimeStep_(0),
    heatOfReactionIonisation_(readScalar(propsDict_.lookup("heatOfReactionIonisation"))),
    heatOfReactionDissociation_(readScalar(propsDict_.lookup("heatOfReactionDissociation"))),
    relax_(true),
    allowSplitting_(false),
    writeRatesToTerminal_(false),
    volume_(0.0),
    numberDensities_(2, 0.0)
{

}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

dissociationIonisationTypeISameSpecies::~dissociationIonisationTypeISameSpecies()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void dissociationIonisationTypeISameSpecies::initialConfiguration()
{
    setProperties();
}

void dissociationIonisationTypeISameSpecies::setProperties()
{
    // reading in reactants

    const List<word> reactantMolecules (propsDict_.lookup("reactantMolecules"));

    if(reactantMolecules.size() != 2)
    {
        FatalErrorIn("dissociationIonisationTypeISameSpecies::setProperties()")
            << "There should be two reactants, instead of " 
            << reactantMolecules.size() << nl 
            << exit(FatalError);
    }
    
    if(reactantMolecules[0] != reactantMolecules[1])
    {
        FatalErrorIn("dissociationIonisationTypeISameSpecies::setProperties()")
            << "Both reactant species must be the same, they are currently " 
	    << reactantMolecules[0] << " and " << reactantMolecules[1] << nl 
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
            FatalErrorIn("dissociationIonisationTypeISameSpecies::setProperties()")
                << "Cannot find type id: " << reactantMolecules[r] << nl 
                << exit(FatalError);
        }

        // check that reactants are 'MOLECULES' (not 'ATOMS') 

        const label& rDof = cloud_.constProps(reactantIds_[r]).rotationalDegreesOfFreedom();
    
        if(rDof < 1)
        {
            FatalErrorIn("dissociationIonisationTypeISameSpecies::setProperties()")
                << "Reactant must be a molecule (not an atom): " << reactantMolecules[r] 
                << nl 
                << exit(FatalError);
        }
        
        const label& vDof = cloud_.constProps(reactantIds_[r]).vibrationalDegreesOfFreedom();
        
        if(vDof > 1)
        {
            FatalErrorIn("dissociationIonisationTypeISameSpecies::setProperties()")
                << "Reactions are currently only implemented for monatomic and diatomic species"
                << " This is a polyatomic:" << reactantMolecules[r] 
                << nl 
                << exit(FatalError);
        }
    }
    
    // reading in ionisation products

    const List<word> productMoleculesIonisation (propsDict_.lookup("productsOfIonisedMolecule"));

    if(productMoleculesIonisation.size() != 2)
    {
        FatalErrorIn("dissociationIonisationTypeISameSpecies::setProperties()")
            << "There should be two products, instead of " 
            << productMoleculesIonisation.size() << nl 
            << exit(FatalError);
    }
    
    productIdsIonisation_.setSize(productMoleculesIonisation.size(), -1);

    forAll(productIdsIonisation_, r)
    {
        productIdsIonisation_[r] = findIndex(cloud_.typeIdList(), productMoleculesIonisation[r]);

        // check that products belong to the typeIdList (constant/dsmcProperties)
        if(productIdsIonisation_[r] == -1)
        {
            FatalErrorIn("dissociationIonisationTypeISameSpecies::setProperties()")
                << "Cannot find type id: " << productMoleculesIonisation[r] << nl 
                << exit(FatalError);
        }
    }
    
    const scalar& rDof1 = cloud_.constProps(productIdsIonisation_[0]).rotationalDegreesOfFreedom();

    if(rDof1 < VSMALL)
    {
        FatalErrorIn("dissociationIonisationTypeISameSpecies::setProperties()")
            << "First ionisation product must be a molecular ion: " << productMoleculesIonisation[0] 
            << nl 
            << exit(FatalError);
    }
    
    const label& charge = cloud_.constProps(productIdsIonisation_[1]).charge();

    if(charge != -1)
    {
        FatalErrorIn("dissociationIonisationTypeISameSpecies::setProperties()")
            << "Second ionisation product must be an electron: " << productMoleculesIonisation[1] 
            << nl 
            << exit(FatalError);
    }

    // reading in dissociation products

    const List<word> productMoleculesDissociation (propsDict_.lookup("productsOfDissociatedMolecule"));

    if(productMoleculesDissociation.size() != 2)
    {
        FatalErrorIn("dissociationIonisationTypeISameSpecies::setProperties()")
            << "There should be two products, instead of " 
            << productMoleculesDissociation.size() << nl 
            << exit(FatalError);
    }
    
    productIdsDissociation_.setSize(productMoleculesDissociation.size(), -1);

    forAll(productIdsDissociation_, r)
    {
        productIdsDissociation_[r] = findIndex(cloud_.typeIdList(), productMoleculesDissociation[r]);

        // check that products belong to the typeIdList (constant/dsmcProperties)
        if(productIdsDissociation_[r] == -1)
        {
            FatalErrorIn("dissociationIonisationTypeISameSpecies::setProperties()")
                << "Cannot find type id: " << productMoleculesDissociation[r] << nl 
                << exit(FatalError);
        }

        // check that products are 'ATOMS' (not 'MOLECULES') 

        const scalar& rDof = cloud_.constProps(productIdsDissociation_[r]).rotationalDegreesOfFreedom();
    
        if(rDof > 1)
        {
            FatalErrorIn("dissociationIonisationTypeISameSpecies::setProperties()")
                << "Reactant must be an atom (not a molecule): " << productMoleculesDissociation[r] 
                << nl 
                << exit(FatalError);
        }
    }
}

bool dissociationIonisationTypeISameSpecies::tryReactMolecules(const label& typeIdP, const label& typeIdQ)
{
    label reactantPId = findIndex(reactantIds_, typeIdP);
    label reactantQId = findIndex(reactantIds_, typeIdQ);

    if(reactantPId == reactantQId)
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

void dissociationIonisationTypeISameSpecies::reaction
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


void dissociationIonisationTypeISameSpecies::reaction
(
    dsmcParcel& p,
    dsmcParcel& q
)
{
    const label typeIdP = p.typeId();
    const label typeIdQ = q.typeId();
    
    if(typeIdP == typeIdQ && typeIdP == reactantIds_[0]) // same species and desired species to measure rate for
    { 
        relax_ = true;
        
        scalar totalReactionProbability = 0.0;
        scalarList reactionProbabilities(2, 0.0);
        
        vector UP = p.U();
        vector UQ = q.U();
        
        scalar ERotP = p.ERot();
        scalar ERotQ = q.ERot();
        const scalar EVibP = p.vibLevel()[0]*cloud_.constProps(typeIdP).thetaV()[0]*physicoChemical::k.value();
        const scalar EVibQ = q.vibLevel()[0]*cloud_.constProps(typeIdQ).thetaV()[0]*physicoChemical::k.value();
        const scalar EEleP = cloud_.constProps(typeIdP).electronicEnergyList()[p.ELevel()];
        const scalar EEleQ = cloud_.constProps(typeIdQ).electronicEnergyList()[q.ELevel()];

        const scalar mP = cloud_.constProps(typeIdP).mass();
        const scalar mQ = cloud_.constProps(typeIdQ).mass();
        
        const scalar thetaVP = cloud_.constProps(typeIdP).thetaV()[0];
        const scalar thetaVQ = cloud_.constProps(typeIdQ).thetaV()[0];
        
        const scalar thetaDQ = cloud_.constProps(typeIdQ).thetaD()[0];
        
        const label idP = cloud_.constProps(typeIdP).charDissQuantumLevel()[0];
        const label idQ = cloud_.constProps(typeIdQ).charDissQuantumLevel()[0];

        const scalar ZrefQ = cloud_.constProps(typeIdQ).Zref()[0];

        const scalar refTempZvQ = cloud_.constProps(typeIdQ).TrefZv()[0];
        
        const scalar rotationalDofQ = cloud_.constProps(typeIdQ).rotationalDegreesOfFreedom();
        
        const List<label>& gListP = cloud_.constProps(typeIdP).degeneracyList();
        const List<scalar>& EElistP = cloud_.constProps(typeIdP).electronicEnergyList();
        
        const label jMaxQ = cloud_.constProps(typeIdQ).numberOfElectronicLevels();
        const List<label>& gListQ = cloud_.constProps(typeIdQ).degeneracyList();
        const List<scalar>& EElistQ = cloud_.constProps(typeIdQ).electronicEnergyList();

        const scalar mR = mP*mQ/(mP + mQ);
        const scalar cRsqr = magSqr(UP - UQ);
        scalar translationalEnergy = 0.5*mR*cRsqr;

        const scalar heatOfReactionIonisationJoules = heatOfReactionIonisation_*physicoChemical::k.value();
        const scalar heatOfReactionDissociationJoules = heatOfReactionDissociation_*physicoChemical::k.value();

        const scalar omegaPQ =
            0.5
            *(
                    cloud_.constProps(typeIdP).omega()
                + cloud_.constProps(typeIdQ).omega()
            );
       
        scalar ChiB = 2.5 - omegaPQ;
        
        //Test for P reactions first
        
        bool dissocReactionP = false;
        bool ionisationReactionP = false;
                
        //2 reactions possible
        // 1. Dissociation of P
        // 2. Ionisation of P

        scalar EcPP = 0.0;
        
        EcPP = translationalEnergy + EVibP;
        label imaxP = EcPP/(physicoChemical::k.value()*thetaVP);
        
        if(imaxP-idP > 0)
        {
            //Dissociation can occur
            totalReactionProbability += 1.0;
            reactionProbabilities[0] = 1.0;
        }
        
        scalar ionisationEnergy = cloud_.constProps(typeIdP).ionisationTemperature()*physicoChemical::k.value();
        
        // calculate if an ionisation of species P is possible
        EcPP = translationalEnergy + EEleP;

        if((EcPP - ionisationEnergy) > VSMALL)
        {
            totalReactionProbability += 1.0;
            reactionProbabilities[1] = 1.0;
        }
        
        //Decide if a reaction is to occur
        
        if(totalReactionProbability > cloud_.rndGen().sample01<scalar>())
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
                    
                    if(cumulativeProbability > cloud_.rndGen().sample01<scalar>())
                    {
                        //Current reaction is to occur
                        
                        if(i == 0)
                        {
                            //Ionisation is to occur
                            dissocReactionP = true;
                            break;
                        }
                        if(i == 1)
                        {
                            //Dissociation reaction is to occur
                            ionisationReactionP = true;
                            break;
                        }
                    }
                }
            }
        }
        
        bool dissocReactionQ = false;
        bool ionisationReactionQ = false;
        
        scalar totalReactionProbabilityQ = 0.0;
        scalarList reactionProbabilitiesQ(2, 0.0);
        
        //2 reactions possible
        // 1. Dissociation of Q
        // 2. Ionisation of Q
        
        scalar EcPQ = translationalEnergy + EVibQ;
        
        if(dissocReactionP)
        {
            EcPQ += heatOfReactionDissociationJoules;
        }
        
        if(ionisationReactionP)
        {
            EcPQ += heatOfReactionIonisationJoules;
        }
        
        label imaxQ = EcPQ/(physicoChemical::k.value()*thetaVQ);
        
        if(imaxQ-idQ > 0)
        {
            //Dissociation can occur
            totalReactionProbabilityQ += 1.0;
            reactionProbabilitiesQ[0] = 1.0;
        }
        
        scalar ionisationEnergyQ = cloud_.constProps(typeIdQ).ionisationTemperature()*physicoChemical::k.value();
        
        // calculate if an ionisation of species Q is possible
        EcPQ = translationalEnergy + EEleQ;
        
        if(dissocReactionP)
        {
            EcPQ += heatOfReactionDissociationJoules;
        }
        
        if(ionisationReactionP)
        {
            EcPQ += heatOfReactionIonisationJoules;
        }

        if((EcPQ - ionisationEnergyQ) > VSMALL)
        {
            totalReactionProbabilityQ += 1.0;
            reactionProbabilitiesQ[1] = 1.0;
        }
        
        //Decide if a reaction is to occur
        
        if(totalReactionProbabilityQ > cloud_.rndGen().sample01<scalar>())
        {
            //A chemical reaction is to occur, choose which one
            
            scalarList normalisedProbabilities(reactionProbabilitiesQ.size(), 0.0);
            scalar cumulativeProbability = 0.0;
            
            normalisedProbabilities = reactionProbabilitiesQ/totalReactionProbabilityQ;
            
            forAll(normalisedProbabilities, i)
            {                
                //If current reaction can't occur, don't check for it
                if(normalisedProbabilities[i] > VSMALL)
                {
                    cumulativeProbability += normalisedProbabilities[i];
                    
                    if(cumulativeProbability > cloud_.rndGen().sample01<scalar>())
                    {
                        //Current reaction is to occur
                        
                        if(i == 0)
                        {
                            //Ionisation is to occur
                            dissocReactionQ = true;
                            break;
                        }
                        if(i == 1)
                        {
                            //Dissociation reaction is to occur
                            ionisationReactionQ = true;
                            break;
                        }
                    }
                }
            }
        }
        
        
        //dissociation of P only
        if(dissocReactionP && !dissocReactionQ &&!ionisationReactionP &&!ionisationReactionQ)
        {
//             nReactionsPerTimeStep_++;
            nTotDissociationReactions_++;
            nDissociationReactionsPerTimeStep_++;
            
            if(allowSplitting_)
            {
                relax_ = false;
                
                translationalEnergy = translationalEnergy + heatOfReactionDissociationJoules + EVibP;
            
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

                const label& typeId1 = productIdsDissociation_[0];
                const label& typeId2 = productIdsDissociation_[1];
                
                //Mass of Product one and two
                scalar mP1 = cloud_.constProps(typeId1).mass();
                scalar mP2 = cloud_.constProps(typeId2).mass();
                
                scalar mRatoms = mP1*mP2/(mP1 + mP2);
                
                translationalEnergy = ERotP + EEleP;
                
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

                // Q remains NON-DISSOCIATED.
                q.U() = UQ;
                q.ERot() = ERotQ;
                q.vibLevel()[0] = vibLevelQ;
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
        
        //dissociation of Q only
        if(dissocReactionQ && !dissocReactionP &&!ionisationReactionP &&!ionisationReactionQ)
        {
//             nReactionsPerTimeStep_++;
            nTotDissociationReactions_++;
            nDissociationReactionsPerTimeStep_++;
            
            if(allowSplitting_)
            {
                relax_ = false;
                
                translationalEnergy = translationalEnergy + heatOfReactionDissociationJoules + EVibQ;
            
                translationalEnergy += EEleP;
                    
                label ELevelP = cloud_.postCollisionElectronicEnergyLevel
                                (
                                    translationalEnergy,
                                    jMaxQ,
                                    omegaPQ,
                                    EElistQ,
                                    gListQ
                                );
                                
                translationalEnergy -= EElistP[ELevelP];
                
                translationalEnergy += EVibP;
                
                label iMax = (translationalEnergy / (physicoChemical::k.value()*thetaVQ));
                
                label vibLevelP = cloud_.postCollisionVibrationalEnergyLevel
                                (
                                        true,
                                        p.vibLevel()[0],
                                        iMax,
                                        thetaVQ,
                                        thetaDQ,
                                        refTempZvQ,
                                        omegaPQ,
                                        ZrefQ,
                                        translationalEnergy
                                    );
                                
                translationalEnergy -= vibLevelP*thetaVQ*physicoChemical::k.value();
                                
                translationalEnergy += ERotP;
                
                ERotP = translationalEnergy*cloud_.postCollisionRotationalEnergy(rotationalDofQ,ChiB);
                        
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

                const label& typeId1 = productIdsDissociation_[0];
                const label& typeId2 = productIdsDissociation_[1];
                
                //Mass of Product one and two
                scalar mP1 = cloud_.constProps(typeId1).mass();
                scalar mP2 = cloud_.constProps(typeId2).mass();
                
                scalar mRatoms = mP1*mP2/(mP1 + mP2);
                
                translationalEnergy = ERotQ + EEleQ;
                
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

                // P remains NON-DISSOCIATED.
                p.U() = UP;
                p.ERot() = ERotP;
                p.vibLevel()[0] = vibLevelP;
                p.ELevel() = ELevelP;

                // Molecule Q will dissociate
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
        
        //dissociation of both P and Q 
        if(dissocReactionP && dissocReactionQ)
        {
//             nReactionsPerTimeStep_++;
            nTotDissociationReactions_++;
            nTotDissociationReactions_++;
            nDissociationReactionsPerTimeStep_++;
            nDissociationReactionsPerTimeStep_++;
            
            if(allowSplitting_)
            {
                relax_ = false;
                
                translationalEnergy = translationalEnergy + 2.0*heatOfReactionDissociationJoules + EVibQ + EVibP;
           
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

                const label& typeId1 = productIdsDissociation_[0];
                const label& typeId2 = productIdsDissociation_[1];
                
                //Mass of Product one and two
                scalar mP1 = cloud_.constProps(typeId1).mass();
                scalar mP2 = cloud_.constProps(typeId2).mass();
                
                scalar mRatoms = mP1*mP2/(mP1 + mP2);
                
                scalar translationalEnergy1 = ERotP + EEleP;
                
                scalar cRatoms1 = sqrt(2.0*translationalEnergy1/mRatoms);
                
                scalar translationalEnergy2 = ERotQ + EEleQ;
                
                scalar cRatoms2 = sqrt(2.0*translationalEnergy2/mRatoms);

                // Variable Hard Sphere collision part
                scalar cosTheta2 = 2.0*cloud_.rndGen().sample01<scalar>() - 1.0;
            
                scalar sinTheta2 = sqrt(1.0 - cosTheta2*cosTheta2);
            
                scalar phi2 = twoPi*cloud_.rndGen().sample01<scalar>();
                
                scalar cosTheta3 = 2.0*cloud_.rndGen().sample01<scalar>() - 1.0;
            
                scalar sinTheta3 = sqrt(1.0 - cosTheta3*cosTheta3);
            
                scalar phi3 = twoPi*cloud_.rndGen().sample01<scalar>();
            
                vector postCollisionRelU1 = cRatoms1
                *vector
                    (
                        cosTheta2,
                        sinTheta2*cos(phi2),
                        sinTheta2*sin(phi2)
                    );
                    
                vector postCollisionRelU2 = cRatoms2
                *vector
                    (
                        cosTheta3,
                        sinTheta3*cos(phi3),
                        sinTheta3*sin(phi3)
                    );

                    
                vector uP1 = UP + postCollisionRelU1*mP2/(mP1 + mP2);
                vector uP2 = UP - postCollisionRelU1*mP1/(mP1 + mP2);
                
                vector uQ1 = UQ + postCollisionRelU2*mP2/(mP1 + mP2);
                vector uQ2 = UQ - postCollisionRelU2*mP1/(mP1 + mP2);

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

                // Molecule Q will dissociate
                vector position2 = q.position();
                
                label cell2 = -1;
                label tetFace2 = -1;
                label tetPt2 = -1;

                mesh_.findCellFacePt
                (
                    position2,
                    cell2,
                    tetFace2,
                    tetPt2
                );

                q.typeId() = typeId1;
                q.U() = uQ1;
                q.vibLevel().setSize(0,0);
                q.ERot() = 0.0;
                q.ELevel() = 0;
                
                label classificationQ = q.classification();
                scalar RWF2 = q.RWF();
                labelList vibLevel2;
                
                // insert new product 2
                cloud_.addNewParcel
                (
                    position2,
                    uQ2,
                    RWF2,
                    0.0,
                    0,
                    cell2,
                    tetFace2,
                    tetPt2,
                    typeId2,
                    -1,
                    classificationQ,
                    vibLevel2
                );
            }
        }

        //ionisation of P only
        if(ionisationReactionP &&!ionisationReactionQ &&!dissocReactionP &&!dissocReactionQ)
        {
//             nReactionsPerTimeStep_++;
            nTotIonisationReactions_++;
            nIonisationReactionsPerTimeStep_++;
            
            if(allowSplitting_)
            {
                relax_ = false;
                
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

                const label& typeId1 = productIdsIonisation_[0];
                const label& typeId2 = productIdsIonisation_[1];
                
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
                q.ERot() = ERotQ;
                q.vibLevel()[0] = vibLevelQ;
                q.ELevel() = ELevelQ;

                // Molecule P will ionise
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
                
                forAll(p.vibLevel(), i)
                {
                    p.vibLevel()[i] = 0;
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
        
       // ionisation of Q only
        if(ionisationReactionQ &&!ionisationReactionP &&!dissocReactionP &&!dissocReactionQ)
        {
//             nReactionsPerTimeStep_++;
            nTotIonisationReactions_++;
            nIonisationReactionsPerTimeStep_++;
            
            if(allowSplitting_)
            {
                relax_ = false;
                
                translationalEnergy = translationalEnergy + heatOfReactionIonisationJoules + EEleQ;
                
                translationalEnergy += EEleP;
                    
                label ELevelP = cloud_.postCollisionElectronicEnergyLevel
                                (
                                    translationalEnergy,
                                    jMaxQ,
                                    omegaPQ,
                                    EElistQ,
                                    gListQ
                                );
                                
                translationalEnergy -= EElistP[ELevelP];
                
                translationalEnergy += EVibP;
                
                label iMax = (translationalEnergy / (physicoChemical::k.value()*thetaVQ));
                
                label vibLevelP = cloud_.postCollisionVibrationalEnergyLevel
                                (
                                        true,
                                        p.vibLevel()[0],
                                        iMax,
                                        thetaVQ,
                                        thetaDQ,
                                        refTempZvQ,
                                        omegaPQ,
                                        ZrefQ,
                                        translationalEnergy
                                    );
                                
                translationalEnergy -= vibLevelP*thetaVQ*physicoChemical::k.value();
                                
                translationalEnergy += ERotP;
                
                ERotP = translationalEnergy*cloud_.postCollisionRotationalEnergy(rotationalDofQ,ChiB);
                        
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

                const label& typeId1 = productIdsIonisation_[0];
                const label& typeId2 = productIdsIonisation_[1];
                
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
                p.ERot() = ERotP;
                p.vibLevel()[0] = vibLevelP;
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
                
                forAll(q.vibLevel(), i)
                {
                    q.vibLevel()[i] = 0;
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
        
        //ionisation of both P and Q 
        if(ionisationReactionP &&ionisationReactionQ)
        {
            nTotIonisationReactions_++;
            nTotIonisationReactions_++;
            nIonisationReactionsPerTimeStep_++;
            nIonisationReactionsPerTimeStep_++;
            
            if(allowSplitting_)
            {
                relax_ = false;
                
                translationalEnergy = translationalEnergy + 2.0*heatOfReactionIonisationJoules + EEleP + EEleQ;

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

                const label& typeId1 = productIdsIonisation_[0];
                const label& typeId2 = productIdsIonisation_[1];
                
                //Mass of Product one and two
                scalar mP1 = cloud_.constProps(typeId1).mass();
                scalar mP2 = cloud_.constProps(typeId2).mass();
                
                scalar mRatoms = mP1*mP2/(mP1 + mP2);
                
                scalar translationalEnergy1 = ERotP + EVibP;
                scalar translationalEnergy2 = ERotQ + EVibQ;
                
                scalar cRatoms1 = sqrt(2.0*translationalEnergy1/mRatoms);
                scalar cRatoms2 = sqrt(2.0*translationalEnergy2/mRatoms);

                // Variable Hard Sphere collision part
                scalar cosTheta2 = 2.0*cloud_.rndGen().sample01<scalar>() - 1.0;
            
                scalar sinTheta2 = sqrt(1.0 - cosTheta2*cosTheta2);
            
                scalar phi2 = twoPi*cloud_.rndGen().sample01<scalar>();
                
                scalar cosTheta3 = 2.0*cloud_.rndGen().sample01<scalar>() - 1.0;
            
                scalar sinTheta3 = sqrt(1.0 - cosTheta3*cosTheta3);
            
                scalar phi3 = twoPi*cloud_.rndGen().sample01<scalar>();
            
                vector postCollisionRelU1 = cRatoms1
                *vector
                    (
                        cosTheta2,
                        sinTheta2*cos(phi2),
                        sinTheta2*sin(phi2)
                    );
                    
                vector postCollisionRelU2 = cRatoms2
                *vector
                    (
                        cosTheta3,
                        sinTheta3*cos(phi3),
                        sinTheta3*sin(phi3)
                    );


                vector uP1 = UP + postCollisionRelU1*mP2/(mP1 + mP2);
                vector uP2 = UP - postCollisionRelU1*mP1/(mP1 + mP2);
                
                vector uQ1 = UQ + postCollisionRelU2*mP2/(mP1 + mP2);
                vector uQ2 = UQ - postCollisionRelU2*mP1/(mP1 + mP2);

                // Molecule P will ionise
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
                
                forAll(p.vibLevel(), i)
                {
                    p.vibLevel()[i] = 0;
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
                
                // Molecule Q will ionise
                vector position2 = q.position();
                
                label cell2 = -1;
                label tetFace2 = -1;
                label tetPt2 = -1;

                mesh_.findCellFacePt
                (
                    position2,
                    cell2,
                    tetFace2,
                    tetPt2
                );
                
                q.typeId() = typeId1;
                q.U() = uQ1;
                
                forAll(q.vibLevel(), i)
                {
                    q.vibLevel()[i] = 0;
                }

                q.ERot() = 0.0;
                q.ELevel() = 0;
                
                label classificationQ = q.classification();
                scalar RWF2 = q.RWF();
                labelList vibLevel2;
                
                // insert new product 2
                cloud_.addNewParcel
                (
                    position2,
                    uQ2,
                    RWF2,
                    0.0,
                    0,
                    cell2,
                    tetFace2,
                    tetPt2,
                    typeId2,
                    -1,
                    classificationQ,
                    vibLevel2
                );
            }
        }
        
        //dissociation P and ionisation of Q 
        if(dissocReactionP &&ionisationReactionQ &&!ionisationReactionP &&!dissocReactionQ)
        {
            nTotIonisationReactions_++;
            nTotDissociationReactions_++;
            nIonisationReactionsPerTimeStep_++;
            nDissociationReactionsPerTimeStep_++;
            
            if(allowSplitting_)
            {
                relax_ = false;
                
                translationalEnergy = translationalEnergy + heatOfReactionIonisationJoules + heatOfReactionDissociationJoules + EVibP + EEleQ;

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

                const label& typeId1 = productIdsDissociation_[0];
                const label& typeId2 = productIdsDissociation_[1];
                
                const label& typeId3 = productIdsIonisation_[0];
                const label& typeId4 = productIdsIonisation_[1];
                
                //Mass of Product one and two
                scalar mP1 = cloud_.constProps(typeId1).mass();
                scalar mP2 = cloud_.constProps(typeId2).mass();
                
                scalar mP3 = cloud_.constProps(typeId3).mass();
                scalar mP4 = cloud_.constProps(typeId4).mass();
                
                scalar mRatoms1 = mP1*mP2/(mP1 + mP2);
                scalar mRatoms2 = mP3*mP4/(mP3 + mP4);
                
                scalar translationalEnergy1 = ERotP + EEleP;
                scalar translationalEnergy2 = ERotQ + EVibQ;
                
                scalar cRatoms1 = sqrt(2.0*translationalEnergy1/mRatoms1);
                scalar cRatoms2 = sqrt(2.0*translationalEnergy2/mRatoms2);

                // Variable Hard Sphere collision part
                scalar cosTheta2 = 2.0*cloud_.rndGen().sample01<scalar>() - 1.0;
            
                scalar sinTheta2 = sqrt(1.0 - cosTheta2*cosTheta2);
            
                scalar phi2 = twoPi*cloud_.rndGen().sample01<scalar>();
                
                scalar cosTheta3 = 2.0*cloud_.rndGen().sample01<scalar>() - 1.0;
            
                scalar sinTheta3 = sqrt(1.0 - cosTheta3*cosTheta3);
            
                scalar phi3 = twoPi*cloud_.rndGen().sample01<scalar>();
            
                vector postCollisionRelU1 = cRatoms1
                *vector
                    (
                        cosTheta2,
                        sinTheta2*cos(phi2),
                        sinTheta2*sin(phi2)
                    );
                    
                vector postCollisionRelU2 = cRatoms2
                *vector
                    (
                        cosTheta3,
                        sinTheta3*cos(phi3),
                        sinTheta3*sin(phi3)
                    );


                vector uP1 = UP + postCollisionRelU1*mP2/(mP1 + mP2);
                vector uP2 = UP - postCollisionRelU1*mP1/(mP1 + mP2);
                
                vector uQ1 = UQ + postCollisionRelU2*mP4/(mP3 + mP4);
                vector uQ2 = UQ - postCollisionRelU2*mP3/(mP3 + mP4);

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
                
                // Molecule Q will ionise
                vector position2 = q.position();
                
                label cell2 = -1;
                label tetFace2 = -1;
                label tetPt2 = -1;

                mesh_.findCellFacePt
                (
                    position2,
                    cell2,
                    tetFace2,
                    tetPt2
                );
                
                q.typeId() = typeId3;
                q.U() = uQ1;
                
                forAll(q.vibLevel(), i)
                {
                    q.vibLevel()[i] = 0;
                }

                q.ERot() = 0.0;
                q.ELevel() = 0;
                
                label classificationQ = q.classification();
                scalar RWF2 = q.RWF();
                labelList vibLevel2;
                
                // insert new product 2
                cloud_.addNewParcel
                (
                    position2,
                    uQ2,
                    RWF2,
                    0.0,
                    0,
                    cell2,
                    tetFace2,
                    tetPt2,
                    typeId4,
                    -1,
                    classificationQ,
                    vibLevel2
                );
            }
        }
        
        //ionisation of P and dissociation of Q 
        if(dissocReactionQ &&ionisationReactionP &&!ionisationReactionP &&!dissocReactionQ)
        {
            nTotIonisationReactions_++;
            nTotDissociationReactions_++;
            nIonisationReactionsPerTimeStep_++;
            nDissociationReactionsPerTimeStep_++;
            
            if(allowSplitting_)
            {
                relax_ = false;
                
                translationalEnergy = translationalEnergy + heatOfReactionIonisationJoules + heatOfReactionDissociationJoules + EEleP + EVibQ;

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

                const label& typeId1 = productIdsIonisation_[0];
                const label& typeId2 = productIdsIonisation_[1];
                
                const label& typeId3 = productIdsDissociation_[0];
                const label& typeId4 = productIdsDissociation_[1];
                
                //Mass of Product one and two
                scalar mP1 = cloud_.constProps(typeId1).mass();
                scalar mP2 = cloud_.constProps(typeId2).mass();
                
                scalar mP3 = cloud_.constProps(typeId3).mass();
                scalar mP4 = cloud_.constProps(typeId4).mass();
                
                scalar mRatoms1 = mP1*mP2/(mP1 + mP2);
                scalar mRatoms2 = mP3*mP4/(mP3 + mP4);
                
                scalar translationalEnergy1 = ERotP + EEleP;
                scalar translationalEnergy2 = ERotQ + EVibQ;
                
                scalar cRatoms1 = sqrt(2.0*translationalEnergy1/mRatoms1);
                scalar cRatoms2 = sqrt(2.0*translationalEnergy2/mRatoms2);

                // Variable Hard Sphere collision part
                scalar cosTheta2 = 2.0*cloud_.rndGen().sample01<scalar>() - 1.0;
            
                scalar sinTheta2 = sqrt(1.0 - cosTheta2*cosTheta2);
            
                scalar phi2 = twoPi*cloud_.rndGen().sample01<scalar>();
                
                scalar cosTheta3 = 2.0*cloud_.rndGen().sample01<scalar>() - 1.0;
            
                scalar sinTheta3 = sqrt(1.0 - cosTheta3*cosTheta3);
            
                scalar phi3 = twoPi*cloud_.rndGen().sample01<scalar>();
            
                vector postCollisionRelU1 = cRatoms1
                *vector
                    (
                        cosTheta2,
                        sinTheta2*cos(phi2),
                        sinTheta2*sin(phi2)
                    );
                    
                vector postCollisionRelU2 = cRatoms2
                *vector
                    (
                        cosTheta3,
                        sinTheta3*cos(phi3),
                        sinTheta3*sin(phi3)
                    );


                vector uP1 = UP + postCollisionRelU1*mP2/(mP1 + mP2);
                vector uP2 = UP - postCollisionRelU1*mP1/(mP1 + mP2);
                
                vector uQ1 = UQ + postCollisionRelU2*mP4/(mP3 + mP4);
                vector uQ2 = UQ - postCollisionRelU2*mP3/(mP3 + mP4);

                // Molecule P will ionise
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
                
                forAll(p.vibLevel(), i)
                {
                    p.vibLevel()[i] = 0;
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
                
                // Molecule Q will dissociate
                vector position2 = q.position();
                
                label cell2 = -1;
                label tetFace2 = -1;
                label tetPt2 = -1;

                mesh_.findCellFacePt
                (
                    position2,
                    cell2,
                    tetFace2,
                    tetPt2
                );
                
                q.typeId() = typeId3;
                q.U() = uQ1;
                q.vibLevel().setSize(0,0);
                q.ERot() = 0.0;
                q.ELevel() = 0;
                
                label classificationQ = q.classification();
                scalar RWF2 = q.RWF();
                labelList vibLevel2;
                
                // insert new product 2
                cloud_.addNewParcel
                (
                    position2,
                    uQ2,
                    RWF2,
                    0.0,
                    0,
                    cell2,
                    tetFace2,
                    tetPt2,
                    typeId4,
                    -1,
                    classificationQ,
                    vibLevel2
                );
            }
        }
    }
}

void  dissociationIonisationTypeISameSpecies::outputResults(const label& counterIndex)
{
    if(writeRatesToTerminal_ == true)
    {
        // measure density 

        const List< DynamicList<dsmcParcel*> >& cellOccupancy
            = cloud_.cellOccupancy();
            
        volume_ = 0.0;

        label molsReactants = 0;

        forAll(cellOccupancy, c)
        {
            const List<dsmcParcel*>& parcelsInCell = cellOccupancy[c];

            forAll(parcelsInCell, pIC)
            {
                dsmcParcel* p = parcelsInCell[pIC];

                if(findIndex(reactantIds_, p->typeId()) != -1)
                {
                    molsReactants++;
                }
            }

            volume_ += mesh_.cellVolumes()[c];
        }
        
        scalar volume = volume_;
        label nTotReactions = nTotReactions_;

        //- Parallel communication
        if(Pstream::parRun())
        {
            reduce(molsReactants, sumOp<label>());
            reduce(volume, sumOp<scalar>());
            reduce(nTotReactions, sumOp<label>());
        }

        numberDensities_[0] = (molsReactants*cloud().nParticle())/volume;
        numberDensities_[1] = (molsReactants*cloud().nParticle())/volume; 

        const scalar& deltaT = mesh_.time().deltaT().value();

        word reactantMolA = cloud_.typeIdList()[reactantIds_[0]];
        word reactantMolB = cloud_.typeIdList()[reactantIds_[1]];

        word productMolA = cloud_.typeIdList()[productIdsDissociation_[0]];
        word productMolB = cloud_.typeIdList()[productIdsDissociation_[1]];
        
        word productMolC = cloud_.typeIdList()[productIdsIonisation_[0]];
        word productMolD = cloud_.typeIdList()[productIdsIonisation_[1]];

        if((numberDensities_[0] > 0.0) && (numberDensities_[1] > 0.0))
        {
            scalar reactionRateIonisation = 0.0;
            scalar reactionRateDissociation = 0.0;
            
            reactionRateIonisation =
            (
                nTotIonisationReactions_
                * cloud_.nParticle()
            )/(counterIndex*deltaT*numberDensities_[0]* numberDensities_[1]*volume);
               
            Info<< "Ionisation type I reaction " 
                <<  reactantMolA << " + " << reactantMolB
                <<  " --> " 
                << productMolC << " + " << productMolD << " + " << reactantMolB 
                << ", reaction rate = " << reactionRateIonisation
                << endl;
                
            reactionRateDissociation =
            (
                nTotDissociationReactions_
                * cloud_.nParticle()
            )/(counterIndex*deltaT*numberDensities_[0]* numberDensities_[1]*volume);
               
            Info<< "Dissociation type I reaction " 
                <<  reactantMolA << " + " << reactantMolB
                <<  " --> " 
                << productMolA << " + " << productMolB << " + " << reactantMolB 
                << ", reaction rate = " << reactionRateDissociation
                << endl;
        }
    }
    else
    {
        label nTotDissociationReactions = nTotDissociationReactions_;   
        label nTotIonisationReactions = nTotIonisationReactions_;
        label nDissociationReactionsPerTimeStep = nDissociationReactionsPerTimeStep_;
        label nIonisationReactionsPerTimeStep = nIonisationReactionsPerTimeStep_;
        
        if(Pstream::parRun())
        {
            reduce(nTotDissociationReactions, sumOp<label>());
            reduce(nTotIonisationReactions, sumOp<label>());
            reduce(nDissociationReactionsPerTimeStep, sumOp<label>());
            reduce(nIonisationReactionsPerTimeStep, sumOp<label>());
        }
       
       if(nTotDissociationReactions > VSMALL)
       {
            word reactantMolA = cloud_.typeIdList()[reactantIds_[0]];
            word reactantMolB = cloud_.typeIdList()[reactantIds_[1]];

            word productMolA = cloud_.typeIdList()[productIdsDissociation_[0]];
            word productMolB = cloud_.typeIdList()[productIdsDissociation_[1]];
           
            Info<< "Dissociation type I reaction " 
                <<  reactantMolA << " + " << reactantMolB
                <<  " --> " 
                << productMolA << " + " << productMolB << " + " << reactantMolB 
                << " is active, nReactions this time step = " << nDissociationReactionsPerTimeStep << endl;
        } 
        
        if(nTotIonisationReactions > VSMALL)
        {
            word reactantMolA = cloud_.typeIdList()[reactantIds_[0]];
            word reactantMolB = cloud_.typeIdList()[reactantIds_[1]];

            word productMolA = cloud_.typeIdList()[productIdsIonisation_[0]];
            word productMolB = cloud_.typeIdList()[productIdsIonisation_[1]];
        
            Info<< "Ionisation type I reaction " 
                <<  reactantMolA << " + " << reactantMolB
                <<  " --> " 
                << productMolA << " + " << productMolB << " + " << reactantMolB 
                << " is active, nReactions this time step = " << nIonisationReactionsPerTimeStep << endl;
        } 
    }

//     nReactionsPerTimeStep_ = 0.0;
    nDissociationReactionsPerTimeStep_ = 0;
    nIonisationReactionsPerTimeStep_ = 0;
}


const bool& dissociationIonisationTypeISameSpecies::relax() const
{
    return relax_;
}

}
// End namespace Foam

// ************************************************************************* //
