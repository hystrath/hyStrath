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

#include "dissociationIonisationTypeIDissimilarSpecies.H"
#include "addToRunTimeSelectionTable.H"
#include "fvc.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(dissociationIonisationTypeIDissimilarSpecies, 0);

addToRunTimeSelectionTable(dsmcReaction, dissociationIonisationTypeIDissimilarSpecies, dictionary);



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
dissociationIonisationTypeIDissimilarSpecies::dissociationIonisationTypeIDissimilarSpecies
(
    Time& t,
    dsmcCloud& cloud,
    const dictionary& dict
)
:
    dsmcReaction(t, cloud, dict),
    propsDict_(dict.subDict(typeName + "Properties")),
    reactantIds_(),
    dissociationProducts_(),
    ionisationProducts_(),
    reactionName_(propsDict_.lookup("reactionName")),
    heatOfReactionDissociationAB_(readScalar(propsDict_.lookup("heatOfReactionDissociationAB"))),
    heatOfReactionIonisationAB_(readScalar(propsDict_.lookup("heatOfReactionIonisationAB"))),
    heatOfReactionDissociationCD_(readScalar(propsDict_.lookup("heatOfReactionDissociationCD"))),
    heatOfReactionIonisationCD_(readScalar(propsDict_.lookup("heatOfReactionIonisationCD"))),
    nTotABDissociationReactions_(0),
    nABDissociationReactionsPerTimeStep_(0),
    nTotABIonisationReactions_(0),
    nABIonisationReactionsPerTimeStep_(0),
    nTotCDDissociationReactions_(0),
    nCDDissociationReactionsPerTimeStep_(0),
    nTotCDIonisationReactions_(0),
    nCDIonisationReactionsPerTimeStep_(0),
    relax_(true),
    allowSplitting_(false),
    writeRatesToTerminal_(false),
    volume_(0.0),
    numberDensities_(2, 0.0)
{

}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

dissociationIonisationTypeIDissimilarSpecies::~dissociationIonisationTypeIDissimilarSpecies()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void dissociationIonisationTypeIDissimilarSpecies::initialConfiguration()
{
    setProperties();
}

void dissociationIonisationTypeIDissimilarSpecies::setProperties()
{
    // reading in reactants

    const List<word> reactantMolecules (propsDict_.lookup("reactantMolecules"));

    if(reactantMolecules.size() != 2)
    {
        FatalErrorIn("dissociationIonisationTypeIDissimilarSpecies::setProperties()")
            << "There should be two or more reactants, instead of " 
            << reactantMolecules.size() << nl 
            << exit(FatalError);
    }
    
    if(reactantMolecules[0] == reactantMolecules[1])
    {
        FatalErrorIn("dissociationIonisationTypeIDissimilarSpecies::setProperties()")
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
            FatalErrorIn("dissociationIonisationTypeIDissimilarSpecies::setProperties()")
                << "Cannot find type id: " << reactantMolecules[r] << nl 
                << exit(FatalError);
        }

        // check that reactants are 'MOLECULES' (not 'ATOMS') 

        const label& rDof = cloud_.constProps(reactantIds_[r]).rotationalDegreesOfFreedom();
    
        if(rDof < 1)
        {
            FatalErrorIn("dissociationIonisationTypeIDissimilarSpecies::setProperties()")
                << "Reactant must be a molecule (not an atom): " << reactantMolecules[r] 
                << nl 
                << exit(FatalError);
        }
        
        const label& vDof = cloud_.constProps(reactantIds_[r]).vibrationalDegreesOfFreedom();
    
        if(vDof > 1)
        {
            FatalErrorIn("dissociationIonisationTypeIDissimilarSpecies::setProperties()")
                << "Reactions are currently only implemented for monatomic and diatomic species"
                << " This is a polyatomic:" << reactantMolecules[r] 
                << nl 
                << exit(FatalError);
        }
    }

    // reading in dissociation products

    List< List<word> > dissociationProducts (propsDict_.lookup("productsOfDissociatedMolecule"));

    if(dissociationProducts.size() != reactantIds_.size())
    {
        FatalErrorIn("dissociationIonisationTypeIDissimilarSpecies::setProperties()")
            << "number of reactant molecules to be dissociated = " << reactantIds_.size()
            << " is not the same as the number of products = " << dissociationProducts.size()
            << exit(FatalError);
    }
    

    dissociationProducts_.setSize(dissociationProducts.size());

    forAll(dissociationProducts_, r)
    {
        const List<word>& productsForDiss = dissociationProducts[r];

        if(productsForDiss.size() != 2)
        {
            FatalErrorIn("dissociationIonisationTypeIDissimilarSpecies::setProperties()")
                << "There should be two products (for the dissociating molecule "
                << reactantMolecules[r] << "), instead of " 
                << productsForDiss.size() << nl 
                << exit(FatalError);
        }
    
        dissociationProducts_[r].setSize(productsForDiss.size(), -1);
    
        forAll(dissociationProducts_[r], p)
        {
            dissociationProducts_[r][p] = findIndex(cloud_.typeIdList(), productsForDiss[p]);
        
            if(dissociationProducts_[r][p] == -1)
            {
                FatalErrorIn("dissociationIonisationTypeIDissimilarSpecies::setProperties()")
                    << "Cannot find type id: " << productsForDiss[p] << nl 
                    << exit(FatalError);
            }
        }
    }
    
    // reading in ionisation products

    List< List<word> > ionisationProducts (propsDict_.lookup("productsOfIonisedMolecule"));

    if(ionisationProducts.size() != reactantIds_.size())
    {
        FatalErrorIn("dissociationIonisationTypeIDissimilarSpecies::setProperties()")
            << "number of reactant molecules to be ionised = " << reactantIds_.size()
            << " is not the same as the number of products = " << ionisationProducts.size()
            << exit(FatalError);
    }
    

    ionisationProducts_.setSize(ionisationProducts.size());

    forAll(ionisationProducts_, r)
    {
        const List<word>& productsForIon = ionisationProducts[r];

        if(productsForIon.size() != 2)
        {
            FatalErrorIn("dissociationIonisationTypeIDissimilarSpecies::setProperties()")
                << "There should be two products (for the ionising molecule "
                << reactantMolecules[r] << "), instead of " 
                << productsForIon.size() << nl 
                << exit(FatalError);
        }
    
        ionisationProducts_[r].setSize(productsForIon.size(), -1);
    
        forAll(ionisationProducts_[r], p)
        {
            ionisationProducts_[r][p] = findIndex(cloud_.typeIdList(), productsForIon[p]);
        
            if(ionisationProducts_[r][p] == -1)
            {
                FatalErrorIn("dissociationIonisationTypeIDissimilarSpecies::setProperties()")
                    << "Cannot find type id: " << productsForIon[p] << nl 
                    << exit(FatalError);
            }
        }
        
        //check ionisation product 2 is an electron
        forAll(ionisationProducts_[r], p)
        {
            if(p == 1)
            {
                const label& charge = cloud_.constProps(ionisationProducts_[r][p]).charge();
                
                if(charge != -1)
                {
                    FatalErrorIn("dissociationIonisationTypeIDissimilarSpecies::setProperties()")
                        << "Second ionisation product must be an electron: " << productsForIon[p] << nl 
                        << exit(FatalError);
                }
            }
        }
    }
}



bool dissociationIonisationTypeIDissimilarSpecies::tryReactMolecules(const label& typeIdP, const label& typeIdQ)
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

void dissociationIonisationTypeIDissimilarSpecies::reaction
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
//£££££££££££££££££££££££££££££££££££££££££££££££££££££££££££££££££££££££££££££££££££££££££££
void dissociationIonisationTypeIDissimilarSpecies::reaction
(
    dsmcParcel& p,
    dsmcParcel& q
)

{
    label typeIdP = p.typeId();
    label typeIdQ = q.typeId();

    if(typeIdP == reactantIds_[0] && typeIdQ == reactantIds_[1])
    {
        scalar totalReactionProbability = 0.0;
        scalarList reactionProbabilities(4, 0.0);
        
        relax_ = true;
        
        vector UP = p.U();
        vector UQ = q.U();
        scalar ERotP = p.ERot();
        scalar ERotQ = q.ERot();
        scalar EVibP = p.vibLevel()[0]*cloud_.constProps(typeIdP).thetaV()[0]*physicoChemical::k.value();
        scalar EVibQ = q.vibLevel()[0]*cloud_.constProps(typeIdQ).thetaV()[0]*physicoChemical::k.value();
        scalar EEleP = cloud_.constProps(typeIdP).electronicEnergyList()[p.ELevel()];
        scalar EEleQ = cloud_.constProps(typeIdQ).electronicEnergyList()[q.ELevel()];

        scalar mP = cloud_.constProps(typeIdP).mass();
        scalar mQ = cloud_.constProps(typeIdQ).mass();

        scalar mR = mP*mQ/(mP + mQ);
        scalar cRsqr = magSqr(UP - UQ);
        scalar translationalEnergy = 0.5*mR*cRsqr;
        
        scalar thetaVP = cloud_.constProps(typeIdP).thetaV()[0];
        scalar thetaVQ = cloud_.constProps(typeIdQ).thetaV()[0];
        
        scalar thetaDP = cloud_.constProps(typeIdP).thetaD()[0];
        scalar thetaDQ = cloud_.constProps(typeIdQ).thetaD()[0];
        
        scalar ZrefP = cloud_.constProps(typeIdP).Zref()[0];
        scalar ZrefQ = cloud_.constProps(typeIdQ).Zref()[0];
        
        scalar refTempZvP = cloud_.constProps(typeIdP).TrefZv()[0];
        scalar refTempZvQ = cloud_.constProps(typeIdQ).TrefZv()[0];
        
        label jMaxP = cloud_.constProps(typeIdP).numberOfElectronicLevels();
        List<label> gListP = cloud_.constProps(typeIdP).degeneracyList();
        List<scalar> EElistP = cloud_.constProps(typeIdP).electronicEnergyList();
        
        label jMaxQ = cloud_.constProps(typeIdQ).numberOfElectronicLevels();
        List<label> gListQ = cloud_.constProps(typeIdQ).degeneracyList();
        List<scalar> EElistQ = cloud_.constProps(typeIdQ).electronicEnergyList();
        
        scalar rotationalDofP = cloud_.constProps(typeIdP).rotationalDegreesOfFreedom();
        scalar rotationalDofQ = cloud_.constProps(typeIdQ).rotationalDegreesOfFreedom();

        scalar omegaPQ =
        0.5
        *(
            cloud_.constProps(typeIdP).omega()
            + cloud_.constProps(typeIdQ).omega()
        );
        
        scalar ChiB = 2.5 - omegaPQ;

        bool dissocReactionP = false;
        bool ionisationReactionP = false;
        bool dissocReactionQ = false;
        bool ionisationReactionQ = false;
                
        //4 reactions possible
        // 1. Dissociation of P
        // 2. Ionisation of P
        // 3. Dissociation of Q
        // 4. Ionisations of Q

        scalar EcPP = 0.0;
        label idP = cloud_.constProps(typeIdP).charDissQuantumLevel()[0];
        label imaxP = 0;

        // calculate if a dissociation of species P is possible
        EcPP = translationalEnergy + EVibP;

        imaxP = EcPP/(physicoChemical::k.value()*cloud_.constProps(typeIdP).thetaV()[0]);
        
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
            //Ionisation can occur
            totalReactionProbability += 1.0;
            reactionProbabilities[1] = 1.0;
        }

        scalar EcQP = 0.0;
        label idQ = cloud_.constProps(typeIdQ).charDissQuantumLevel()[0];
        label imaxQ = 0;

        // calculate if a dissociation of species Q is possible
        EcQP = translationalEnergy + EVibQ;
        
        imaxQ = EcQP/(physicoChemical::k.value()*cloud_.constProps(typeIdQ).thetaV()[0]);

        if(imaxQ-idQ > 0) 
        {
            //Dissociation can occur
            totalReactionProbability += 1.0;
            reactionProbabilities[2] = 1.0;
        }
        
        ionisationEnergy = cloud_.constProps(typeIdQ).ionisationTemperature()*physicoChemical::k.value();
        
        // calculate if an ionisation of species Q is possible
        EcQP = translationalEnergy + EEleQ;

        if((EcQP - ionisationEnergy) > VSMALL)
        {
            //Ionisation can occur
            totalReactionProbability += 1.0;
            reactionProbabilities[3] = 1.0;
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
                            dissocReactionQ = true;
                            break;
                        }
                        if(i == 3)
                        {
                            //Ionisation reaction is to occur
                            ionisationReactionQ = true;
                            break;
                        }
                    }
                }
            }
        }
        
        if(dissocReactionP)
        {
//             nReactionsPerTimeStep_++;
            nTotABDissociationReactions_++;
            nABDissociationReactionsPerTimeStep_++;
            
            if(allowSplitting_)
            {
                relax_ = false;
                
                const scalar& heatOfReactionDissociationJoules = heatOfReactionDissociationAB_*physicoChemical::k.value();
                
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

                const label& typeId1 = dissociationProducts_[0][0];
                const label& typeId2 = dissociationProducts_[0][1];
                
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

                // Molecule P will dissociate.
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
        
        if(dissocReactionQ)
        {
//             nReactionsPerTimeStep_++;
            nTotCDDissociationReactions_++;
            nCDDissociationReactionsPerTimeStep_++;
            
            if(allowSplitting_)
            {
                relax_ = false;
                
                const scalar& heatOfReactionDissociationJoules = heatOfReactionDissociationCD_*physicoChemical::k.value();
                
                translationalEnergy = translationalEnergy + heatOfReactionDissociationJoules + EVibQ;
                
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

                UP = Ucm + (postCollisionRelU*mQ/(mP + mQ)); // UP is the single atom
                UQ = Ucm - (postCollisionRelU*mP/(mP + mQ)); // UQ is used as Ucm for atomic split.

                const label& typeId1 = dissociationProducts_[1][0];
                const label& typeId2 = dissociationProducts_[1][1];
                
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
        
        if(ionisationReactionP)
        {
//             nReactionsPerTimeStep_++;
            nTotABIonisationReactions_++;
            nABIonisationReactionsPerTimeStep_++;
            
            if(allowSplitting_)
            {
                relax_ = false;
                
                const scalar& heatOfReactionIonisationJoules = heatOfReactionIonisationAB_*physicoChemical::k.value();
                
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

                const label& typeId1 = ionisationProducts_[0][0];
                const label& typeId2 = ionisationProducts_[0][1];
                
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

                // Molecule P will dissociate.
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
        
        if(ionisationReactionQ)
        {
//             nReactionsPerTimeStep_++;
            nTotCDIonisationReactions_++;
            nCDIonisationReactionsPerTimeStep_++;
            
            if(allowSplitting_)
            {
                relax_ = false;
                
                const scalar& heatOfReactionIonisationJoules = heatOfReactionIonisationCD_*physicoChemical::k.value();
                
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

                UP = Ucm + (postCollisionRelU*mQ/(mP + mQ)); // UP is the single atom
                UQ = Ucm - (postCollisionRelU*mP/(mP + mQ)); // UQ is used as Ucm for atomic split.

                const label& typeId1 = ionisationProducts_[1][0];
                const label& typeId2 = ionisationProducts_[1][1];
                
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
                q.vibLevel().setSize(0,0);
                q.ERot() = 0.0;
                q.ELevel() = 0;
                
                label classificationP = q.classification();
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
                    classificationP,
                    vibLevel
                );
            }
        }
    }

//////////////////////////////////////////////////////////////////////////////////////////////
  
    if(typeIdP == reactantIds_[1] && typeIdQ == reactantIds_[0])
    {
        relax_ = true;
        
        scalar totalReactionProbability = 0.0;
        scalarList reactionProbabilities(4, 0.0);
        
        vector UP = p.U();
        vector UQ = q.U();
        scalar ERotP = p.ERot();
        scalar ERotQ = q.ERot();
        scalar EVibP = p.vibLevel()[0]*cloud_.constProps(typeIdP).thetaV()[0]*physicoChemical::k.value();
        scalar EVibQ = q.vibLevel()[0]*cloud_.constProps(typeIdQ).thetaV()[0]*physicoChemical::k.value();
        scalar EEleP = cloud_.constProps(typeIdP).electronicEnergyList()[p.ELevel()];
        scalar EEleQ = cloud_.constProps(typeIdQ).electronicEnergyList()[q.ELevel()];

        scalar mP = cloud_.constProps(typeIdP).mass();
        scalar mQ = cloud_.constProps(typeIdQ).mass();

        scalar mR = mP*mQ/(mP + mQ);
        scalar cRsqr = magSqr(UP - UQ);
        scalar translationalEnergy = 0.5*mR*cRsqr;
        
        scalar thetaVP = cloud_.constProps(typeIdP).thetaV()[0];
        scalar thetaVQ = cloud_.constProps(typeIdQ).thetaV()[0];
        
        scalar thetaDP = cloud_.constProps(typeIdP).thetaD()[0];
        scalar thetaDQ = cloud_.constProps(typeIdQ).thetaD()[0];
        
        scalar ZrefP = cloud_.constProps(typeIdP).Zref()[0];
        scalar ZrefQ = cloud_.constProps(typeIdQ).Zref()[0];
        
        scalar refTempZvP = cloud_.constProps(typeIdP).TrefZv()[0];
        scalar refTempZvQ = cloud_.constProps(typeIdQ).TrefZv()[0];
        
        scalar rotationalDofP = cloud_.constProps(typeIdP).rotationalDegreesOfFreedom();
        scalar rotationalDofQ = cloud_.constProps(typeIdQ).rotationalDegreesOfFreedom();
        
        label jMaxP = cloud_.constProps(typeIdP).numberOfElectronicLevels();
        List<label> gListP = cloud_.constProps(typeIdP).degeneracyList();
        List<scalar> EElistP = cloud_.constProps(typeIdP).electronicEnergyList();
        
        label jMaxQ = cloud_.constProps(typeIdQ).numberOfElectronicLevels();
        List<label> gListQ = cloud_.constProps(typeIdQ).degeneracyList();
        List<scalar> EElistQ = cloud_.constProps(typeIdQ).electronicEnergyList();

        scalar omegaPQ =
        0.5
        *(
            cloud_.constProps(typeIdP).omega()
            + cloud_.constProps(typeIdQ).omega()
        );
        
        scalar ChiB = 2.5 - omegaPQ;

        bool dissocReactionP = false;
        bool ionisationReactionP = false;
        bool dissocReactionQ = false;
        bool ionisationReactionQ = false;
                
        //4 reactions possible
        // 1. Dissociation of P
        // 2. Ionisation of P
        // 3. Dissociation of Q
        // 4. Ionisations of Q

        scalar EcPQ = 0.0;
        label idP = cloud_.constProps(typeIdP).charDissQuantumLevel()[0];
        label imaxP = 0;

        // calculate if a dissociation of species P is possible
        EcPQ = translationalEnergy + EVibP;

        imaxP = EcPQ/(physicoChemical::k.value()*cloud_.constProps(typeIdP).thetaV()[0]);
        
        if(imaxP-idP > 0)
        {
            //Dissociation can occur
            totalReactionProbability += 1.0;
            reactionProbabilities[0] = 1.0;
        }
        
        scalar ionisationEnergy = cloud_.constProps(typeIdP).ionisationTemperature()*physicoChemical::k.value();
        
        // calculate if an ionisation of species P is possible
        EcPQ = translationalEnergy + EEleP;

        if((EcPQ - ionisationEnergy) > VSMALL)
        {
            totalReactionProbability += 1.0;
            reactionProbabilities[1] = 1.0;
        }

        scalar EcQP = 0.0;
        label idQ = cloud_.constProps(typeIdQ).charDissQuantumLevel()[0];
        label imaxQ = 0;

        // calculate if a dissociation of species Q is possible
        EcQP = translationalEnergy + EVibQ;

        imaxQ = EcQP/(physicoChemical::k.value()*cloud_.constProps(typeIdQ).thetaV()[0]);
        
        if(imaxQ-idQ > 0) 
        {
            //Dissociation can occur
            totalReactionProbability += 1.0;
            reactionProbabilities[2] = 1.0;
        }
        
        ionisationEnergy = cloud_.constProps(typeIdQ).ionisationTemperature()*physicoChemical::k.value();
        
        // calculate if an ionisation of species Q is possible
        EcQP = translationalEnergy + EEleQ;

        if((EcQP - ionisationEnergy) > VSMALL)
        {
            totalReactionProbability += 1.0;
            reactionProbabilities[3] = 1.0;
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
//                 normalisedProbabilities[i] = reactionProbabilities[i]/totalReactionProbability;
                
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
                        if(i == 2)
                        {
                            //Ionisation is to occur
                            dissocReactionQ = true;
                            break;
                        }
                        if(i == 3)
                        {
                            //Dissociation reaction is to occur
                            ionisationReactionQ = true;
                            break;
                        }
                    }
                }
            }
        }
        
        if(dissocReactionP)
        {
//             nReactionsPerTimeStep_++;
            nTotCDDissociationReactions_++;
            nCDDissociationReactionsPerTimeStep_++;
            
            if(allowSplitting_)
            {
                relax_ = false;
                
                const scalar& heatOfReactionDissociationJoules = heatOfReactionDissociationCD_*physicoChemical::k.value();
                
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

                const label typeId1 = dissociationProducts_[1][0];
                const label typeId2 = dissociationProducts_[1][1];
                
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

                // Molecule P will dissociate.
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
        
        if(dissocReactionQ)
        {
//             nReactionsPerTimeStep_++;
            nTotABDissociationReactions_++;
            nABDissociationReactionsPerTimeStep_++;
            
            if(allowSplitting_)
            {
                relax_ = false;
                
                const scalar& heatOfReactionDissociationJoules = heatOfReactionDissociationAB_*physicoChemical::k.value();
                
                translationalEnergy = translationalEnergy + heatOfReactionDissociationJoules + EVibQ;
                
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

                UP = Ucm + (postCollisionRelU*mQ/(mP + mQ)); // UP is the single atom
                UQ = Ucm - (postCollisionRelU*mP/(mP + mQ)); // UQ is used as Ucm for atomic split.

                const label typeId1 = dissociationProducts_[0][0];
                const label typeId2 = dissociationProducts_[0][1];
                
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
        
        if(ionisationReactionP)
        {
//             nReactionsPerTimeStep_++;
            nTotCDIonisationReactions_++;
            nCDIonisationReactionsPerTimeStep_++;
            
            if(allowSplitting_)
            {
                relax_ = false;
                
                const scalar& heatOfReactionIonisationJoules = heatOfReactionIonisationCD_*physicoChemical::k.value();
                
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

                const label& typeId1 = ionisationProducts_[1][0];
                const label& typeId2 = ionisationProducts_[1][1];
                
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

                // Q remains NON-DISSOCIATED.
                q.U() = UQ;
                q.ERot() = ERotQ;
                q.vibLevel()[0] = vibLevelQ;
                q.ELevel() = ELevelQ;

                // Molecule P will dissociate.
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
        
        if(ionisationReactionQ)
        {
//             nReactionsPerTimeStep_++;
            nTotABIonisationReactions_++;
            nABIonisationReactionsPerTimeStep_++;
            
            if(allowSplitting_)
            {
                relax_ = false;
                
                const scalar& heatOfReactionIonisationJoules = heatOfReactionIonisationAB_*physicoChemical::k.value();
                
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

                UP = Ucm + (postCollisionRelU*mQ/(mP + mQ)); // UP is the single atom
                UQ = Ucm - (postCollisionRelU*mP/(mP + mQ)); // UQ is used as Ucm for atomic split.

                const label& typeId1 = ionisationProducts_[0][0];
                const label& typeId2 = ionisationProducts_[0][1];
                
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

                // Q remains NON-DISSOCIATED.
                p.U() = UP;
                p.ERot() = ERotP;
                p.vibLevel()[0] = vibLevelP;
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
                
                q.typeId() = typeId1;
                q.U() = uQ1;
                q.vibLevel().setSize(0,0);
                q.ERot() = 0.0;
                q.ELevel() = 0;
                
                label classificationP = q.classification();
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
                    classificationP,
                    vibLevel
                );
            }
        }
    }
}
//&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
void  dissociationIonisationTypeIDissimilarSpecies::outputResults(const label& counterIndex)
{  
    if(writeRatesToTerminal_ == true)
    {
        // measure density 
        const List< DynamicList<dsmcParcel*> >& cellOccupancy
            = cloud_.cellOccupancy();
            
        volume_ = 0.0;

        List<label> mols;
	      mols.append(0); mols.append(0);

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

            volume_ += mesh_.cellVolumes()[c];
        }
        
        scalar volume = volume_;
        label nTotABDissociationReactions = nTotABDissociationReactions_;
        label nTotCDDissociationReactions = nTotCDDissociationReactions_;
        label nTotABIonisationReactions = nTotABIonisationReactions_;
        label nTotCDIonisationReactions = nTotCDIonisationReactions_;

        //- Parallel communication
        if(Pstream::parRun())
        {
            reduce(volume, sumOp<scalar>());
            reduce(mols[0], sumOp<label>());
            reduce(mols[1], sumOp<label>());
            reduce(nTotABDissociationReactions, sumOp<label>());
            reduce(nTotCDDissociationReactions, sumOp<label>());
            reduce(nTotABIonisationReactions, sumOp<label>());
            reduce(nTotCDIonisationReactions, sumOp<label>());
        }

        numberDensities_[0] = (mols[0]*cloud().nParticle())/volume;
        numberDensities_[1] = (mols[1]*cloud().nParticle())/volume;

        const scalar& deltaT = mesh_.time().deltaT().value();

        word reactantMolA = cloud_.typeIdList()[reactantIds_[0]];
        word reactantMolB = cloud_.typeIdList()[reactantIds_[1]];
        
        word dissociationProductMolA = cloud_.typeIdList()[dissociationProducts_[0][0]];
        word dissociationProductMolB = cloud_.typeIdList()[dissociationProducts_[0][1]];
        word dissociationProductMolC = cloud_.typeIdList()[dissociationProducts_[1][0]];
        word dissociationProductMolD = cloud_.typeIdList()[dissociationProducts_[1][1]];
        
        word ionisationProductMolA = cloud_.typeIdList()[ionisationProducts_[0][0]];
        word ionisationProductMolB = cloud_.typeIdList()[ionisationProducts_[0][1]];
        word ionisationProductMolC = cloud_.typeIdList()[ionisationProducts_[1][0]];
        word ionisationProductMolD = cloud_.typeIdList()[ionisationProducts_[1][1]];

        if((numberDensities_[0] > 0.0) && (numberDensities_[1] > 0.0))
        {
            scalar reactionRate1 = 0.0;
            scalar reactionRate2 = 0.0;
            scalar reactionRate3 = 0.0;
            scalar reactionRate4 = 0.0;

            reactionRate1 =
            (
                nTotABDissociationReactions
                * cloud_.nParticle()
            )/(counterIndex*deltaT*numberDensities_[0]* numberDensities_[1]*volume);
        
            reactionRate2 =
            (
                nTotCDDissociationReactions
                * cloud_.nParticle()
            )/(counterIndex*deltaT*numberDensities_[0]* numberDensities_[1]*volume);
            
            reactionRate3 =
            (
                nTotABIonisationReactions
                * cloud_.nParticle()
            )/(counterIndex*deltaT*numberDensities_[0]* numberDensities_[1]*volume);
            
            reactionRate4 =
            (
                nTotCDIonisationReactions
                * cloud_.nParticle()
            )/(counterIndex*deltaT*numberDensities_[0]* numberDensities_[1]*volume);

            Info << "Dissociation type I reaction " <<  reactantMolA << " + " << reactantMolB << " --> " <<
            dissociationProductMolA << " + " << dissociationProductMolB << " + " << reactantMolB <<
            ", reaction rate = " << reactionRate1  << nl      
            << "Dissociation type I reaction " <<  reactantMolB << " + " << reactantMolA << " --> " <<
            dissociationProductMolC << " + " << dissociationProductMolD << " + " << reactantMolA <<
            ", reaction rate = " << reactionRate2 << nl
            << "Ionisation type I reaction " <<  reactantMolA << " + " << reactantMolB << " --> " <<
            ionisationProductMolA << " + " << ionisationProductMolB << " + " << reactantMolB <<
            ", reaction rate = " << reactionRate3 << nl
            << "Ionisation type I reaction " <<  reactantMolB << " + " << reactantMolA << " --> " <<
            ionisationProductMolC << " + " << ionisationProductMolD << " + " << reactantMolA <<
            ", reaction rate = " << reactionRate4
            << endl;
        }
    }
    else
    {
        label nTotABDissociationReactions = nTotABDissociationReactions_;
        label nTotCDDissociationReactions = nTotCDDissociationReactions_;
        label nTotABIonisationReactions = nTotABIonisationReactions_;
        label nTotCDIonisationReactions = nTotCDIonisationReactions_;
        
        label nABDissociationReactionsPerTimeStep = nABDissociationReactionsPerTimeStep_;
        label nCDDissociationReactionsPerTimeStep = nCDDissociationReactionsPerTimeStep_;
        label nABIonisationReactionsPerTimeStep = nABIonisationReactionsPerTimeStep_;
        label nCDIonisationReactionsPerTimeStep = nCDIonisationReactionsPerTimeStep_;
        
        //- Parallel communication
        if(Pstream::parRun())
        {
            reduce(nTotABDissociationReactions, sumOp<label>());
            reduce(nTotCDDissociationReactions, sumOp<label>());
            reduce(nTotABIonisationReactions, sumOp<label>());
            reduce(nTotCDIonisationReactions, sumOp<label>());
            
            reduce(nABDissociationReactionsPerTimeStep, sumOp<label>());
            reduce(nCDDissociationReactionsPerTimeStep, sumOp<label>());
            reduce(nABIonisationReactionsPerTimeStep, sumOp<label>());
            reduce(nCDIonisationReactionsPerTimeStep, sumOp<label>());
        }
        
        if(nTotABDissociationReactions > VSMALL)
        {
            word reactantMolA = cloud_.typeIdList()[reactantIds_[0]];
            word reactantMolB = cloud_.typeIdList()[reactantIds_[1]];
            
            word dissociationProductMolA = cloud_.typeIdList()[dissociationProducts_[0][0]];
            word dissociationProductMolB = cloud_.typeIdList()[dissociationProducts_[0][1]];
            
            Info << "Dissociation type I reaction " <<  reactantMolA << " + " << reactantMolB << " --> " <<
                dissociationProductMolA << " + " << dissociationProductMolB << " + " << reactantMolB <<
                " is active, nReactions this time step = " << nABDissociationReactionsPerTimeStep << endl;;
        }
        
        if(nTotCDDissociationReactions > VSMALL)
        {
            word reactantMolA = cloud_.typeIdList()[reactantIds_[0]];
            word reactantMolB = cloud_.typeIdList()[reactantIds_[1]];
            
            word dissociationProductMolC = cloud_.typeIdList()[dissociationProducts_[1][0]];
            word dissociationProductMolD = cloud_.typeIdList()[dissociationProducts_[1][1]];
               
            Info << "Dissociation type I reaction " <<  reactantMolB << " + " << reactantMolA << " --> " <<
                dissociationProductMolC << " + " << dissociationProductMolD << " + " << reactantMolA <<
                " is active, nReactions this time step = " << nCDDissociationReactionsPerTimeStep << endl;
        }
        
        if(nTotABIonisationReactions > VSMALL)
        {
            word reactantMolA = cloud_.typeIdList()[reactantIds_[0]];
            word reactantMolB = cloud_.typeIdList()[reactantIds_[1]];
            
            word ionisationProductMolA = cloud_.typeIdList()[ionisationProducts_[0][0]];
            word ionisationProductMolB = cloud_.typeIdList()[ionisationProducts_[0][1]];
            
            Info << "Ionisation type I reaction " <<  reactantMolA << " + " << reactantMolB << " --> " <<
                ionisationProductMolA << " + " << ionisationProductMolB << " + " << reactantMolB <<
                " is active, nReactions this time step = " << nABIonisationReactionsPerTimeStep << endl;
        }
        
        if(nTotCDIonisationReactions > VSMALL)
        {
            word reactantMolA = cloud_.typeIdList()[reactantIds_[0]];
            word reactantMolB = cloud_.typeIdList()[reactantIds_[1]];
            
            word ionisationProductMolC = cloud_.typeIdList()[ionisationProducts_[1][0]];
            word ionisationProductMolD = cloud_.typeIdList()[ionisationProducts_[1][1]];
               
            Info << "Ionisation type I reaction " <<  reactantMolB << " + " << reactantMolA << " --> " <<
                ionisationProductMolC << " + " << ionisationProductMolD << " + " << reactantMolA <<
                " is active, nReactions this time step = " << nCDIonisationReactionsPerTimeStep << endl;
        }

    }

//     nReactionsPerTimeStep_ = 0.0;
    nABDissociationReactionsPerTimeStep_ = 0;
    nCDDissociationReactionsPerTimeStep_ = 0;
    nABIonisationReactionsPerTimeStep_ = 0;
    nCDIonisationReactionsPerTimeStep_ = 0;

}


const bool& dissociationIonisationTypeIDissimilarSpecies::relax() const
{
    return relax_;
}

}
// End namespace Foam

// ************************************************************************* //
