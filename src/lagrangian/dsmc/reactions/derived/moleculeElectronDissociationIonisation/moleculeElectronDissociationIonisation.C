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

#include "moleculeElectronDissociationIonisation.H"
#include "addToRunTimeSelectionTable.H"
#include "fvc.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(moleculeElectronDissociationIonisation, 0);

addToRunTimeSelectionTable(dsmcReaction, moleculeElectronDissociationIonisation, dictionary);



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
moleculeElectronDissociationIonisation::moleculeElectronDissociationIonisation
(
    Time& t,
    dsmcCloud& cloud,
    const dictionary& dict
)
:
    dsmcReaction(t, cloud, dict),
    propsDict_(dict.subDict(typeName + "Properties")),
    reactantIds_(),
    productIdsDiss_(),
    productIdsIon_(),
    reactionName_(propsDict_.lookup("reactionName")),
    heatOfReactionDiss_(readScalar(propsDict_.lookup("heatOfReactionDissociation"))),
    heatOfReactionIon_(readScalar(propsDict_.lookup("heatOfReactionIonisation"))),
    nDissociationReactions_(0),
    nDissociationReactionsPerTimeStep_(0),
    nIonisationReactions_(0),
    nIonisationReactionsPerTimeStep_(0),
    relax_(true),
    allowSplitting_(false),
    writeRatesToTerminal_(false),
    volume_(0.0),
    numberDensities_(2, 0.0)
{

}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

moleculeElectronDissociationIonisation::~moleculeElectronDissociationIonisation()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void moleculeElectronDissociationIonisation::initialConfiguration()
{
    setProperties();
}

void moleculeElectronDissociationIonisation::setProperties()
{
    // reading in reactants

    const List<word> reactantMolecules (propsDict_.lookup("reactants"));

    if(reactantMolecules.size() != 2)
    {
        FatalErrorIn("moleculeElectronDissociationIonisation::setProperties()")
            << "There should be two or more reactants, instead of " 
            << reactantMolecules.size() << nl 
            << exit(FatalError);
    }
    
    if(reactantMolecules[0] == reactantMolecules[1])
    {
        FatalErrorIn("moleculeElectronDissociationIonisation::setProperties()")
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
            FatalErrorIn("moleculeElectronDissociationIonisation::setProperties()")
                << "Cannot find type id: " << reactantMolecules[r] << nl 
                << exit(FatalError);
        }
    }
    
    // check that reactant one is a 'MOLECULE'' 

    const label& rDof1 = cloud_.constProps(reactantIds_[0]).rotationalDegreesOfFreedom();

    if(rDof1 < 1)
    {
        FatalErrorIn("moleculeElectronDissociationIonisation::setProperties()")
            << "First reactant must be a molecule (not an atom or an electron): " << reactantMolecules[0] 
            << nl 
            << exit(FatalError);
    }
    
    const label& vDof = cloud_.constProps(reactantIds_[0]).vibrationalDegreesOfFreedom();

    if(vDof > 1)
    {
        FatalErrorIn("moleculeAtomDissociationIonisation::setProperties()")
            << "Reactions are currently only implemented for monatomic and diatomic species" 
            << " This is a polyatomic:" << reactantMolecules[0] 
            << nl 
            << exit(FatalError);
    }
    
    // check that reactant two is an 'ELECTRON'

    const label& charge = cloud_.constProps(reactantIds_[1]).charge();

    if(charge != -1)
    {
        FatalErrorIn("moleculeElectronDissociationIonisation::setProperties()")
            << "Second reactant must be an electron, not " << reactantMolecules[1] 
            << nl 
            << exit(FatalError);
    }
    
    // reading in dissociation products

    List<word> productMoleculesDissociation (propsDict_.lookup("productsOfDissociatedMolecule"));

    if(productMoleculesDissociation.size() != 2)
    {
        FatalErrorIn("moleculeElectronDissociationIonisation::setProperties()")
            << "Number of dissociation products is " << productMoleculesDissociation.size() <<
            ", should be two."
            << exit(FatalError);
    }
    

    productIdsDiss_.setSize(productMoleculesDissociation.size());

    forAll(productMoleculesDissociation, r)
    {
        if(productIdsDiss_.size() != 2)
        {
            FatalErrorIn("moleculeElectronDissociationIonisation::setProperties()")
                << "There should be two products (for the dissociating molecule "
                << reactantMolecules[r] << "), instead of " 
                << productIdsDiss_.size() << nl 
                << exit(FatalError);
        }
    
        forAll(productIdsDiss_, r)
        {
            productIdsDiss_[r] = findIndex(cloud_.typeIdList(), productMoleculesDissociation[r]);

            // check that reactants belong to the typeIdList (constant/dsmcProperties)
            if(productIdsDiss_[r] == -1)
            {
                FatalErrorIn("moleculeElectronDissociationIonisation::setProperties()")
                    << "Cannot find type id: " << productMoleculesDissociation[r] << nl 
                    << exit(FatalError);
            }
        }
        
        // check that product one is an 'ATOM' (not a 'MOLECULE') 

        const scalar& rDof3 = cloud_.constProps(productIdsDiss_[0]).rotationalDegreesOfFreedom();

        if(rDof3 != 0)
        {
            FatalErrorIn("moleculeElectronDissociationIonisation::setProperties()")
                << "First product must be an atom (not a molecule): " << productMoleculesDissociation[0] 
                << nl 
                << exit(FatalError);
        }

        // check that product two is an 'ATOM' (not a 'MOLECULE') 

        const scalar& rDof4 = cloud_.constProps(productIdsDiss_[1]).rotationalDegreesOfFreedom();

        if(rDof4 != 0)
        {
            FatalErrorIn("moleculeElectronDissociationIonisation::setProperties()")
                << "Second product must be an atom (not a molecule): " << productMoleculesDissociation[1] 
                << nl 
                << exit(FatalError);
        }
    }
    
    // reading in ionisation products

    List<word> productMoleculesIonisation (propsDict_.lookup("productsOfIonisedMolecule"));

    if(productMoleculesIonisation.size() != 2)
    {
        FatalErrorIn("moleculeElectronDissociationIonisation::setProperties()")
            << "Number of dissociation products is " << productMoleculesIonisation.size() <<
            ", should be two."
            << exit(FatalError);
    }
    

    productIdsIon_.setSize(productMoleculesIonisation.size());

    forAll(productMoleculesIonisation, r)
    {
        if(productIdsIon_.size() != 2)
        {
            FatalErrorIn("moleculeElectronDissociationIonisation::setProperties()")
                << "There should be two products (for the dissociating molecule "
                << reactantMolecules[r] << "), instead of " 
                << productIdsIon_.size() << nl 
                << exit(FatalError);
        }
    
        forAll(productIdsIon_, r)
        {
            productIdsIon_[r] = findIndex(cloud_.typeIdList(), productMoleculesIonisation[r]);

            // check that reactants belong to the typeIdList (constant/dsmcProperties)
            if(productIdsIon_[r] == -1)
            {
                FatalErrorIn("moleculeElectronDissociationIonisation::setProperties()")
                    << "Cannot find type id: " << productMoleculesIonisation[r] << nl 
                    << exit(FatalError);
            }
        }
        
        // check that product one is a 'MOLECULE', not an 'ATOM'

        const scalar& rDof5 = cloud_.constProps(productIdsIon_[0]).rotationalDegreesOfFreedom();

        if(rDof5 < 1)
        {
            FatalErrorIn("moleculeElectronDissociationIonisation::setProperties()")
                << "First product must be a molecule (not an atom): " << productMoleculesIonisation[0] 
                << nl 
                << exit(FatalError);
        }

        // check that product two is an 'ELECTRON'

        const label& charge = cloud_.constProps(productIdsIon_[1]).charge();

        if(charge != -1)
        {
            FatalErrorIn("moleculeElectronDissociationIonisation::setProperties()")
                << "Second product must be an electron: " << productMoleculesIonisation[1] 
                << nl 
                << exit(FatalError);
        }
    }
}

bool moleculeElectronDissociationIonisation::tryReactMolecules(const label& typeIdP, const label& typeIdQ)
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

void moleculeElectronDissociationIonisation::reaction
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


void moleculeElectronDissociationIonisation::reaction
(
    dsmcParcel& p,
    dsmcParcel& q
)
{
    label typeIdP = p.typeId();
    label typeIdQ = q.typeId();

    if(typeIdP == reactantIds_[0] && typeIdQ == reactantIds_[1]) // This produces the correct equilibrium rate A2 + X.
    {
        relax_ = true;
        
        scalar totalReactionProbability = 0.0;
        scalarList reactionProbabilities(2, 0.0);
        
        vector UP = p.U();
        vector UQ = q.U();
        scalar ERotP = p.ERot();
        scalar EVibP = p.vibLevel()[0]*cloud_.constProps(typeIdP).thetaV()[0]*physicoChemical::k.value();
        scalar EEleP = cloud_.constProps(typeIdP).electronicEnergyList()[p.ELevel()];

        scalar mP = cloud_.constProps(typeIdP).mass();
        scalar mQ = cloud_.constProps(typeIdQ).mass();

        scalar mR = mP*mQ/(mP + mQ);
        scalar cRsqr = magSqr(UP - UQ);
        scalar translationalEnergy = 0.5*mR*cRsqr;
        
        scalar thetaVP = cloud_.constProps(typeIdP).thetaV()[0];
        List<label> gListP = cloud_.constProps(typeIdP).degeneracyList();
        List<scalar> EElistP = cloud_.constProps(typeIdP).electronicEnergyList();

        scalar heatOfReactionJoulesDiss = heatOfReactionDiss_*physicoChemical::k.value();
        scalar heatOfReactionJoulesIon = heatOfReactionIon_*physicoChemical::k.value();

        
        bool dissocReaction = false;
        bool ionisationReaction = false;
        
        scalar EcPPIon= 0.0;
        scalar ionisationEnergy = cloud_.constProps(typeIdP).ionisationTemperature()*physicoChemical::k.value();
        
        // calculate if an ionisation of species P is possible
        EcPPIon = translationalEnergy + EEleP;
        
        if((EcPPIon - ionisationEnergy) > VSMALL)
        {
            //IONISATION CAN OCCUR
            totalReactionProbability += 1.0;
            reactionProbabilities[0] = 1.0;
        }

        scalar EcPPDiss = 0.0;
        label idP = cloud_.constProps(typeIdP).charDissQuantumLevel()[0];
        label imaxP = 0;
        
        // calculate if a dissociation of species P is possible
        EcPPDiss = translationalEnergy + EVibP;

        imaxP = EcPPDiss/(physicoChemical::k.value()*thetaVP);

        if(imaxP-idP > 0)
        {
            //DISSOCIATION CAN OCCUR
            totalReactionProbability += 1.0;
            reactionProbabilities[1] = 1.0;
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
                            //Ionisation is to occur
                            ionisationReaction = true;
                            break;
                        }
                        if(i == 1)
                        {
                            //Dissociation reaction is to occur
                            dissocReaction = true;
                            break;
                        }
                    }
                }
            }
        }
        
        //Perform a dissociation reaction
        if(dissocReaction)
        {
            nDissociationReactions_++;
            nDissociationReactionsPerTimeStep_++;

            if(allowSplitting_)
            {  
                relax_ = false;
                                
                translationalEnergy = translationalEnergy + heatOfReactionJoulesDiss + EVibP;
                
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
    
                UP = Ucm + postCollisionRelU*mQ/(mP + mQ);
                UQ = Ucm - postCollisionRelU*mP/(mP + mQ); // Q is the NON-DISSOCIATING molecule.

                const label& typeId1 = productIdsDiss_[0];
                const label& typeId2 = productIdsDiss_[1];

                //Mass of Product one and two
                scalar mP1 = cloud_.constProps(typeId1).mass();
                scalar mP2 = cloud_.constProps(typeId2).mass();

                scalar mRatoms = mP1*mP2/(mP1 + mP2);
                
                //center of mass velocity of all particles

                vector UcmAtoms = UP;
                
                label ELevelAtom1 = 0;
                label ELevelAtom2 = 0;
                
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

                vector uP1 = UcmAtoms + postCollisionRelU2*mP2/(mP1 + mP2);
                vector uP2 = UcmAtoms - postCollisionRelU2*mP1/(mP1 + mP2);

                // New electron velocity
                q.U() = UQ;
                q.ELevel() = 0;

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
                p.vibLevel().setSize(0,0);
                p.ERot() = 0.0;
                p.ELevel() = ELevelAtom1;
                
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
                    ELevelAtom2,
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
        
        //Perform an ionisation reaction
        
        if(ionisationReaction)
        {
            nIonisationReactions_++;
            nIonisationReactionsPerTimeStep_++;

            if(allowSplitting_)
            {
                relax_ = false;
                
                label ELevelQ = 0;
                                
                translationalEnergy = translationalEnergy + heatOfReactionJoulesIon + EEleP;
                
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
    
                UP = Ucm + postCollisionRelU*mQ/(mP + mQ);
                UQ = Ucm - postCollisionRelU*mP/(mP + mQ); // Q is the NON-IONISING atom.

                label p1 = 0;
                label p2 = 1;

                const label& typeId1 = productIdsIon_[p1];
                const label& typeId2 = productIdsIon_[p2];
                
                //Mass of Product one and two
                scalar mP1 = cloud_.constProps(typeId1).mass(); // ion
                scalar mP2 = cloud_.constProps(typeId2).mass(); // electron

                scalar mRatoms = mP1*mP2/(mP1 + mP2);
                
                //center of mass velocity of all particles

                vector UcmAtoms = UP;
                
                scalar translationalEnergy2 = ERotP + EVibP;

                scalar cRatoms = sqrt(2.0*translationalEnergy2/mRatoms);

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

                // New electron Q velocity.
                q.U() = UQ;
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
                p.vibLevel().setSize(1,0);
                p.ERot() = 0.0;
                p.ELevel() = 0;
                
                label classificationP = p.classification();
                scalar RWF = p.RWF();
                labelList vibLevel;
		            vibLevel.append(0);
                
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
    }
  
    if(typeIdP == reactantIds_[1] && typeIdQ == reactantIds_[0]) // This produces the correct equilibrium rate X + A2.
    {
        relax_ = true;
        
        scalar totalReactionProbability = 0.0;
        scalarList reactionProbabilities(2, 0.0);
        
        vector UP = p.U();
        vector UQ = q.U();
        scalar ERotQ = q.ERot();
        scalar EVibQ = q.vibLevel()[0]*cloud_.constProps(typeIdQ).thetaV()[0]*physicoChemical::k.value();
        scalar EEleQ = cloud_.constProps(typeIdQ).electronicEnergyList()[q.ELevel()];

        scalar mP = cloud_.constProps(typeIdP).mass();
        scalar mQ = cloud_.constProps(typeIdQ).mass();

        scalar mR = mP*mQ/(mP + mQ);
        scalar cRsqr = magSqr(UP - UQ);
        scalar translationalEnergy = 0.5*mR*cRsqr;
        
        List<label> gListP = cloud_.constProps(typeIdP).degeneracyList();
        List<scalar> EElistP = cloud_.constProps(typeIdP).electronicEnergyList();
        
        List<label> gListQ = cloud_.constProps(typeIdQ).degeneracyList();
        List<scalar> EElistQ = cloud_.constProps(typeIdQ).electronicEnergyList();
        
        scalar heatOfReactionJoulesDiss = heatOfReactionDiss_*physicoChemical::k.value();
        scalar heatOfReactionJoulesIon = heatOfReactionIon_*physicoChemical::k.value();
        
        bool dissocReaction = false;
        bool ionisationReaction = false;
        
        scalar EcPQIon= 0.0;
        scalar ionisationEnergy = cloud_.constProps(typeIdQ).ionisationTemperature()*physicoChemical::k.value();
        
        // calculate if an ionisation of species Q is possible
        EcPQIon = translationalEnergy + EEleQ;
        
        if((EcPQIon - ionisationEnergy) > VSMALL)
        {
            //IONISATION CAN OCCUR
            totalReactionProbability += 1.0;
            reactionProbabilities[0] = 1.0;
        }
        
        scalar thetaVQ = cloud_.constProps(typeIdQ).thetaV()[0];

        scalar EcPQ = 0.0;
        label idQ = cloud_.constProps(typeIdQ).charDissQuantumLevel()[0];
        label imaxQ = 0;

        // calculate if a dissociation of species Q is possible
        EcPQ = translationalEnergy + EVibQ;

        imaxQ = EcPQ/(physicoChemical::k.value()*thetaVQ);	    

        if(imaxQ-idQ > 0)
        {
            //DISSOCIATION CAN OCCUR
            totalReactionProbability += 1.0;
            reactionProbabilities[1] = 1.0;
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
                            //Ionisation is to occur
                            ionisationReaction = true;
                            break;
                        }
                        if(i == 1)
                        {
                            //Dissociation reaction is to occur
                            dissocReaction = true;
                            break;
                        }
                    }
                }
            }
        }
        
        //Perform a dissociation reaction
        if(dissocReaction)
        {
            nDissociationReactions_++;
            nDissociationReactionsPerTimeStep_++;
            
            if(allowSplitting_)
            {
                relax_ = false;
                
                translationalEnergy = translationalEnergy + heatOfReactionJoulesDiss + EVibQ;
                
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
                
                label q1 = 0;
                label q2 = 1;

                const label& typeId1 = productIdsDiss_[q1];
                const label& typeId2 = productIdsDiss_[q2];
                
                //Mass of Product one and two
                scalar mP1 = cloud_.constProps(typeId1).mass();
                scalar mP2 = cloud_.constProps(typeId2).mass();
                
                scalar mRatoms = mP1*mP2/(mP1 + mP2);
                
                label ELevelAtom1 = 0;
                label ELevelAtom2 = 0;
                
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


                vector uP1 = UQ + postCollisionRelU2*mP2/(mP1 + mP2);
                vector uP2 = UQ - postCollisionRelU2*mP1/(mP1 + mP2);

                // P is the electron.
                p.U() = UP;
                p.ELevel() = 0;

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
                q.U() = uP1;
                q.vibLevel().setSize(0,0);
                q.ERot() = 0.0;
                q.ELevel() = ELevelAtom1;
                
                label classificationQ = q.classification();
                scalar RWF = q.RWF();
                labelList vibLevel;
                
                // insert new product 2
                cloud_.addNewParcel
                (
                    position,
                    uP2,
                    RWF,
                    0.0,
                    ELevelAtom2,
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
        
        //Perform an ionisation reaction
        if(ionisationReaction)
        {
            nIonisationReactions_++;
            nIonisationReactionsPerTimeStep_++;
            
            if(allowSplitting_)
            {
                relax_ = false;
                
                label ELevelP = 0;
                                
                translationalEnergy = translationalEnergy + heatOfReactionJoulesIon + EEleQ;
                
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
    
                UP = Ucm + postCollisionRelU*mQ/(mP + mQ);
                UQ = Ucm - postCollisionRelU*mP/(mP + mQ); // Q is the NON-IONISING atom.

                const label& typeId1 = productIdsIon_[0];
                const label& typeId2 = productIdsIon_[1];
                
                //Mass of Product one and two
                scalar mP1 = cloud_.constProps(typeId1).mass(); // ion
                scalar mP2 = cloud_.constProps(typeId2).mass(); // electron

                scalar mRatoms = mP1*mP2/(mP1 + mP2);
                
                //center of mass velocity of all particles

                vector UcmAtoms = UQ;
                
                scalar translationalEnergy2 = ERotQ + EVibQ;

                scalar cRatoms = sqrt(2.0*translationalEnergy2/mRatoms);

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

                // New atom P velocity.
                p.U() = UP;
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
                q.U() = uP1;
                q.vibLevel().setSize(1,0);
                q.ERot() = 0.0;
                q.ELevel() = 0;
                
                label classificationP = q.classification();
                scalar RWF = q.RWF();
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
    }
}

void  moleculeElectronDissociationIonisation::outputResults(const label& counterIndex)
{    
    if(writeRatesToTerminal_ == true)
    {
        const List< DynamicList<dsmcParcel*> >& cellOccupancy
            = cloud_.cellOccupancy();

        List<label> mols;
	mols.append(0); mols.append(0);
        volume_ = 0.0;

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
        label nTotReactionsDissociation = nDissociationReactions_;
        label nTotReactionsIonisation = nIonisationReactions_;

        //- Parallel communication
        if(Pstream::parRun())
        {
            reduce(volume, sumOp<scalar>());
            reduce(mols[0], sumOp<label>());
            reduce(mols[1], sumOp<label>());
            reduce(nTotReactionsDissociation, sumOp<label>());
            reduce(nTotReactionsIonisation, sumOp<label>());
        }

        numberDensities_[0] = (mols[0]*cloud().nParticle())/volume;
        numberDensities_[1] = (mols[1]*cloud().nParticle())/volume;

        const scalar& deltaT = mesh_.time().deltaT().value();

        word reactantMolA = cloud_.typeIdList()[reactantIds_[0]];
        word reactantMolB = cloud_.typeIdList()[reactantIds_[1]];

        word productMolA = cloud_.typeIdList()[productIdsDiss_[0]];
        word productMolB = cloud_.typeIdList()[productIdsDiss_[1]];
        
        word productMolC = cloud_.typeIdList()[productIdsIon_[0]];
        word productMolD = cloud_.typeIdList()[productIdsIon_[1]];

        if((numberDensities_[0] > 0.0) && (numberDensities_[1] > 0.0))
        { 

            scalar reactionRateDissociation = 0.0;

            reactionRateDissociation =
            (
            nTotReactionsDissociation
            * cloud_.nParticle()
            )/(counterIndex*deltaT*numberDensities_[0]* numberDensities_[1]*volume);

            Info<< "Electron dissociation reaction "
                <<  reactantMolA << " + " << reactantMolB 
                <<  " --> "
                << productMolA << " + " << productMolB << " + " << reactantMolB 
                << ", reaction rate = " << reactionRateDissociation
                << endl;
                
            scalar reactionRateIonisation = 0.0;

            reactionRateIonisation =
            (
            nTotReactionsIonisation
            * cloud_.nParticle()
            )/(counterIndex*deltaT*numberDensities_[0]* numberDensities_[1]*volume);

            Info<< "Electron ionisation reaction "
                <<  reactantMolA << " + " << reactantMolB 
                <<  " --> "
                << productMolC << " + " << productMolD << " + " << reactantMolB 
                << ", reaction rate = " << reactionRateIonisation
                << endl;
        }
    }
    else
    {
        label nTotReactionsDissociation = nDissociationReactions_;   
        label nTotReactionsIonisation = nIonisationReactions_;  
        label nDissociationReactionsPerTimeStep = nDissociationReactionsPerTimeStep_;
        label nIonisationReactionsPerTimeStep = nIonisationReactionsPerTimeStep_;
        
        if(Pstream::parRun())
        {
            reduce(nTotReactionsDissociation, sumOp<label>());
            reduce(nTotReactionsIonisation, sumOp<label>());
            
            reduce(nDissociationReactionsPerTimeStep, sumOp<label>());
            reduce(nIonisationReactionsPerTimeStep, sumOp<label>());
        }
    
        if(nTotReactionsDissociation > VSMALL)
        {
                word reactantMolA = cloud_.typeIdList()[reactantIds_[0]];
                word reactantMolB = cloud_.typeIdList()[reactantIds_[1]];

                word productMolA = cloud_.typeIdList()[productIdsDiss_[0]];
                word productMolB = cloud_.typeIdList()[productIdsDiss_[1]];
            
                Info<< "Electron dissociation reaction "
                <<  reactantMolA << " + " << reactantMolB 
                <<  " --> "
                << productMolA << " + " << productMolB << " + " << reactantMolB  
                    << " is active, nReactions this time step = " << nDissociationReactionsPerTimeStep << endl;
        }  
        
        if(nTotReactionsIonisation > VSMALL)
        {
                word reactantMolA = cloud_.typeIdList()[reactantIds_[0]];
                word reactantMolB = cloud_.typeIdList()[reactantIds_[1]];

                word productMolA = cloud_.typeIdList()[productIdsIon_[0]];
                word productMolB = cloud_.typeIdList()[productIdsIon_[1]];
            
                Info<< "Electron ionisation reaction "
                <<  reactantMolA << " + " << reactantMolB 
                <<  " --> "
                << productMolA << " + " << productMolB << " + " << reactantMolB  
                    << " is active, nReactions this time step = " << nIonisationReactionsPerTimeStep << endl;
        } 
    }

    nReactionsPerTimeStep_ = 0.0;
    nDissociationReactionsPerTimeStep_ = 0.0;
    nIonisationReactionsPerTimeStep_ = 0.0;

}


const bool& moleculeElectronDissociationIonisation::relax() const
{
    return relax_;
}

}
// End namespace Foam

// ************************************************************************* //
