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

#include "dissociationTypeISameSpecies.H"
#include "addToRunTimeSelectionTable.H"
#include "fvc.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(dissociationTypeISameSpecies, 0);

addToRunTimeSelectionTable(dsmcReaction, dissociationTypeISameSpecies, dictionary);



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
dissociationTypeISameSpecies::dissociationTypeISameSpecies
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
    nDissociationReactions_(0),
    nIonisationReactions_(0),
    heatOfReactionDiss_(readScalar(propsDict_.lookup("heatOfReactionDissociation"))),
    heatOfReactionIon_(readScalar(propsDict_.lookup("heatOfReactionIonisation"))),
    relax_(true),
    allowSplitting_(false),
    writeRatesToTerminal_(false),
    volume_(0.0),
    numberDensities_(2, 0.0)
{

}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

dissociationTypeISameSpecies::~dissociationTypeISameSpecies()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void dissociationTypeISameSpecies::initialConfiguration()
{
    setProperties();
}

void dissociationTypeISameSpecies::setProperties()
{
    // reading in reactants

    const List<word> reactantMolecules (propsDict_.lookup("reactantMolecules"));

    if(reactantMolecules.size() != 2)
    {
        FatalErrorIn("dissociationTypeISameSpecies::setProperties()")
            << "There should be two reactants, instead of " 
            << reactantMolecules.size() << nl 
            << exit(FatalError);
    }
    
    if(reactantMolecules[0] != reactantMolecules[1])
    {
        FatalErrorIn("dissociationTypeISameSpecies::setProperties()")
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
            FatalErrorIn("dissociationTypeISameSpecies::setProperties()")
                << "Cannot find type id: " << reactantMolecules[r] << nl 
                << exit(FatalError);
        }

        // check that reactants are 'MOLECULES' (not 'ATOMS') 

        const scalar& rDof = cloud_.constProps(reactantIds_[r]).rotationalDegreesOfFreedom();
    
        if(rDof < 1)
        {
            FatalErrorIn("dissociationTypeISameSpecies::setProperties()")
                << "Reactant must be a molecule (not an atom): " << reactantMolecules[r] 
                << nl 
                << exit(FatalError);
        }
    }

    // reading in products

    const List<word> productMoleculesDiss (propsDict_.lookup("productsOfDissociatedMolecule"));

    if(productMoleculesDiss.size() != 2)
    {
        FatalErrorIn("dissociationTypeISameSpecies::setProperties()")
            << "There should be two products, instead of " 
            << productMoleculesDiss.size() << nl 
            << exit(FatalError);
    }
    
    productIdsDiss_.setSize(productMoleculesDiss.size(), -1);

    forAll(productIdsDiss_, r)
    {
        productIdsDiss_[r] = findIndex(cloud_.typeIdList(), productMoleculesDiss[r]);

        // check that products belong to the typeIdList (constant/dsmcProperties)
        if(productIdsDiss_[r] == -1)
        {
            FatalErrorIn("dissociationTypeISameSpecies::setProperties()")
                << "Cannot find type id: " << productMoleculesDiss[r] << nl 
                << exit(FatalError);
        }

        // check that products are 'ATOMS' (not 'MOLECULES') 

        const scalar& rDof = cloud_.constProps(productIdsDiss_[r]).rotationalDegreesOfFreedom();
    
        if(rDof > 1)
        {
            FatalErrorIn("dissociationTypeISameSpecies::setProperties()")
                << "Reactant must be an atom (not a molecule): " << productMoleculesDiss[r] 
                << nl 
                << exit(FatalError);
        }
    }
    
    const List<word> productMoleculesIon (propsDict_.lookup("productsOfIonisedMolecule"));

    if(productMoleculesIon.size() != 2)
    {
        FatalErrorIn("dissociationTypeISameSpecies::setProperties()")
            << "There should be two products, instead of " 
            << productMoleculesIon.size() << nl 
            << exit(FatalError);
    }
    
    productIdsIon_.setSize(productMoleculesIon.size(), -1);

    forAll(productIdsIon_, r)
    {
        productIdsIon_[r] = findIndex(cloud_.typeIdList(), productMoleculesIon[r]);

        // check that products belong to the typeIdList (constant/dsmcProperties)
        if(productIdsIon_[r] == -1)
        {
            FatalErrorIn("dissociationTypeISameSpecies::setProperties()")
                << "Cannot find type id: " << productMoleculesIon[r] << nl 
                << exit(FatalError);
        }
    }
    
    // check that the first product is an ION

    const scalar& rDof = cloud_.constProps(productIdsIon_[0]).rotationalDegreesOfFreedom();

    if(rDof < 1)
    {
        FatalErrorIn("dissociationTypeISameSpecies::setProperties()")
            << "Reactant must be a molecular ion: " << productMoleculesIon[0] 
            << nl 
            << exit(FatalError);
    }
    
    // check that the second product is an ELECTRON
    
    const scalar& mass = cloud_.constProps(productIdsIon_[1]).mass();
    
    if(mass > 1e-24)
    {
        FatalErrorIn("dissociationTypeISameSpecies::setProperties()")
            << "Reactant must be an electron: " << productMoleculesIon[1] 
            << nl 
            << exit(FatalError);
    }
}

bool dissociationTypeISameSpecies::tryReactMolecules(const label& typeIdP, const label& typeIdQ)
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

void dissociationTypeISameSpecies::reaction
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


void dissociationTypeISameSpecies::reaction
(
    dsmcParcel& p,
    dsmcParcel& q
)
{
    label typeIdP = p.typeId();
    label typeIdQ = q.typeId();
    
    relax_ = true;    
    
    if(typeIdP == typeIdQ && typeIdP == reactantIds_[0]) // same species and desired species to measure rate for
    {
        scalar totalReactionProbability = 0.0;
        scalarList reactionProbabilities(2, 0.0);
        
        vector UP = p.U();
        vector UQ = q.U();
        
        scalar ERotP = p.ERot();
//         scalar ERotQ = q.ERot();
        scalar EVibP = p.vibLevel()*cloud_.constProps(typeIdP).thetaV()*physicoChemical::k.value();
        scalar EVibQ = q.vibLevel()*cloud_.constProps(typeIdQ).thetaV()*physicoChemical::k.value();
        scalar EEleP = cloud_.constProps(typeIdP).electronicEnergyList()[p.ELevel()];
        scalar EEleQ = cloud_.constProps(typeIdQ).electronicEnergyList()[q.ELevel()];

        scalar mP = cloud_.constProps(typeIdP).mass();
        scalar mQ = cloud_.constProps(typeIdQ).mass();
        
        scalar thetaVQ = cloud_.constProps(typeIdQ).thetaV();
//         scalar thetaDQ = cloud_.constProps(typeIdQ).thetaD();
        
        label jMaxP = cloud_.constProps(typeIdP).numberOfElectronicLevels();
        label jMaxQ = cloud_.constProps(typeIdQ).numberOfElectronicLevels();
        List<label> gListP = cloud_.constProps(typeIdP).degeneracyList();
        List<label> gListQ = cloud_.constProps(typeIdQ).degeneracyList();
        List<scalar> EElistP = cloud_.constProps(typeIdP).electronicEnergyList();
        List<scalar> EElistQ = cloud_.constProps(typeIdQ).electronicEnergyList();
        
//         scalar ZrefQ = cloud_.constProps(typeIdQ).Zref();

//         scalar refTempZvQ = cloud_.constProps(typeIdQ).TrefZv();
        
//         scalar rotationalDofP = cloud_.constProps(typeIdP).rotationalDegreesOfFreedom();
//         scalar rotationalDofQ = cloud_.constProps(typeIdQ).rotationalDegreesOfFreedom();

        scalar mR = mP*mQ/(mP + mQ);
        scalar cRsqr = magSqr(UP - UQ);
        scalar translationalEnergy = 0.5*mR*cRsqr;

        scalar heatOfReactionJoulesDiss = heatOfReactionDiss_*physicoChemical::k.value();
        scalar heatOfReactionJoulesIon = heatOfReactionDiss_*physicoChemical::k.value();

        scalar omegaPQ =
            0.5
            *(
                    cloud_.constProps(typeIdP).omega()
                + cloud_.constProps(typeIdQ).omega()
            );
       
//         scalar ChiB = 2.5 - omegaPQ;
        
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
        label idP = cloud_.constProps(typeIdP).thetaD()/cloud_.constProps(typeIdP).thetaV();
        label imaxP = 0;
        
        // calculate if a dissociation of species P is possible
        EcPPDiss = translationalEnergy + EVibP;

        imaxP = EcPPDiss/(physicoChemical::k.value()*cloud_.constProps(typeIdP).thetaV());
        
        if((imaxP-idP) > 0)
        {
            //DISSOCIATION CAN OCCUR
            totalReactionProbability += 1.0;
            reactionProbabilities[1] = 1.0;
        }

        //Decide if a reaction is to occur
        
        if(totalReactionProbability > cloud_.rndGen().scalar01())
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
                    
                    if(cumulativeProbability > cloud_.rndGen().scalar01())
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
            
            if(allowSplitting_)
            {
                relax_ = false;
                
                //Particle P is dissociating, so no vibrational redistribution to it
                
                label ELevelP = 0;
                label ELevelQ = 0;
                
                translationalEnergy = translationalEnergy + heatOfReactionJoulesDiss + EVibP + EEleP;

                if(jMaxP > 1)
                {
                    translationalEnergy += EEleP;
                    
                    ELevelP = cloud_.postCollisionElectronicEnergyLevel
                                    (
                                        translationalEnergy,
                                        jMaxP,
                                        omegaPQ,
                                        EElistP,
                                        gListP
                                    );
                                    
                    translationalEnergy -= EElistP[ELevelP];
                }
                
                if(jMaxQ > 1)
                {
                    
                    translationalEnergy += EEleQ;
                    
                    ELevelQ = cloud_.postCollisionElectronicEnergyLevel
                                    (
                                        translationalEnergy,
                                        jMaxQ,
                                        omegaPQ,
                                        EElistQ,
                                        gListQ
                                    );
                                    
                    translationalEnergy -= EElistQ[ELevelQ];
                }
                
                scalar EcQ = translationalEnergy + EVibQ;
                label iMaxQ = (EcQ /(physicoChemical::k.value()*thetaVQ));
                

                label vibLevelQ = cloud_.postCollisionVibrationalEnergyLevel
                            (
                                true,
                                q.vibLevel(),
                                iMaxQ,
                                cloud_.constProps(typeIdQ).thetaV(),
                                cloud_.constProps(typeIdQ).thetaD(),
                                cloud_.constProps(typeIdQ).TrefZv(),
                                omegaPQ,
                                cloud_.constProps(typeIdQ).Zref(),
                                EcQ
                            );
                            
                translationalEnergy -= (vibLevelQ*cloud_.constProps(typeIdQ).thetaV()*physicoChemical::k.value());

                scalar rotationalDofQ = cloud_.constProps(typeIdQ).rotationalDegreesOfFreedom();
                
                scalar energyRatio = cloud_.postCollisionRotationalEnergy(rotationalDofQ,(2.5 - omegaPQ));

                scalar ERotQ = energyRatio*translationalEnergy;
        
                translationalEnergy -= ERotQ;
                
                scalar relVelNonDissoMol = sqrt((2.0*translationalEnergy)/mR);

                // centre of mass velocity of molecules (pre-split)
                vector Ucm = (mP*UP + mQ*UQ)/(mP + mQ);

                // Variable Hard Sphere collision part for collision of molecules
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
        
                UP = Ucm + (postCollisionRelU*mQ/(mP + mQ)); // UP is used as Ucm for atomic split.
                UQ = Ucm - (postCollisionRelU*mP/(mP + mQ)); // Q is the NON-DISSOCIATING molecule.

                // Q remains NON-DISSOCIATED (but internal and translational energy modified).
                q.ELevel() = ELevelQ;
                q.vibLevel() = vibLevelQ;
                q.ERot() = ERotQ;
                q.U() = UQ;

                //split particle P into the atoms
                
                const label& typeId1 = productIdsDiss_[0];
                const label& typeId2 = productIdsDiss_[1];

                //Mass of Product one and two
                scalar mP1 = cloud_.constProps(typeId1).mass();
                scalar mP2 = cloud_.constProps(typeId2).mass();

                scalar mRatoms = mP1*mP2/(mP1 + mP2);

                //centre of mass velocity
                vector UcmAtoms = UP;
                
                //the relative translational energy for the atoms
                //comes from the internal energy of the complex molecule
                
                //first get new electronic energies for the atoms
                
                // collision energy of particle Q = relative translational energy (ERotP) + pre-collision electronic energy
            
                translationalEnergy = ERotP;
                
                label ELevelAtomP = 0;
                label ELevelAtomQ = 0;
                
                label jMaxAtomP = cloud_.constProps(typeId1).numberOfElectronicLevels();
                List<label> gListAtomP = cloud_.constProps(typeId1).degeneracyList();
                List<scalar> EElistAtomP = cloud_.constProps(typeId1).electronicEnergyList();
                
                label jMaxAtomQ = cloud_.constProps(typeId2).numberOfElectronicLevels();
                List<label> gListAtomQ = cloud_.constProps(typeId2).degeneracyList();
                List<scalar> EElistAtomQ = cloud_.constProps(typeId2).electronicEnergyList();
                
                scalar omegaAtoms = 0.5*(cloud_.constProps(typeId1).omega() + cloud_.constProps(typeId2).omega());
                
                if(jMaxAtomP > 1)
                {
//                     translationalEnergy += EElistAtomP[ELevelAtomP];
                    
                    ELevelAtomP = cloud_.postCollisionElectronicEnergyLevel
                                    (
                                        translationalEnergy,
                                        jMaxAtomP,
                                        omegaAtoms,
                                        EElistAtomP,
                                        gListAtomP
                                    );
                                    
                    translationalEnergy -= EElistAtomP[ELevelAtomP];
                }
                
                // Determine if electronic energy level transition can occur for Q and to which level...
                
                // collision energy of particle Q = relative translational energy (ERotP) + pre-collision electronic energy
                
//                 scalar EcAtomQ = translationalEnergy + EElistAtomQ[ELevelAtomQ]; 
                
                // Determine if electronic energy level transition can occur for Q and to which level...
                        
                if(jMaxAtomP > 1)
                {
//                     translationalEnergy += EElistAtomP[ELevelAtomP];
                    
                    ELevelAtomQ = cloud_.postCollisionElectronicEnergyLevel
                                    (
                                        translationalEnergy,
                                        jMaxAtomQ,
                                        omegaAtoms,
                                        EElistAtomQ,
                                        gListAtomQ
                                    );
                                    
                    translationalEnergy -= EElistAtomQ[ELevelAtomQ];
                }
                
                scalar cRatoms = sqrt(2.0*translationalEnergy/mRatoms);
                
                // Variable Hard Sphere collision part for atoms
                cosTheta = 2.0*cloud_.rndGen().scalar01() - 1.0;
            
                sinTheta = sqrt(1.0 - cosTheta*cosTheta);
            
                phi = twoPi*cloud_.rndGen().scalar01();
            
                vector postCollisionRelU2 = cRatoms
                *vector
                    (
                        cosTheta,
                        sinTheta*cos(phi),
                        sinTheta*sin(phi)
                    );

                vector uP1 = UcmAtoms + (postCollisionRelU2*mP2/(mP1 + mP2));
                vector uP2 = UcmAtoms - (postCollisionRelU2*mP1/(mP1 + mP2));

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
                    
                //convert the molecule to an atom
                p.typeId() = typeId1;
                p.U() = uP1;
                p.vibLevel() = 0;
                p.ERot() = 0.0;
                p.ELevel() = ELevelAtomP;
                
                label classificationP = p.classification();
                scalar RWF = p.RWF();

                // insert the second atom
                cloud_.addNewParcel
                (
                    position,
                    uP2,
                    RWF,
                    0.0,
                    0,
                    ELevelAtomQ,
                    cell,
                    tetFace,
                    tetPt,
                    typeId2,
                    0,
                    classificationP
                );
            }
        }
        
        //Perform an ionisation reaction
        if(ionisationReaction)
        {
            nIonisationReactions_++;
        }
    }
}

void  dissociationTypeISameSpecies::outputResults(const label& counterIndex)
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
        label nTotDissociationReactions = nDissociationReactions_;
        label nTotIonisationReactions = nIonisationReactions_;

        //- Parallel communication
        if(Pstream::parRun())
        {
            reduce(molsReactants, sumOp<label>());
            reduce(volume, sumOp<scalar>());
            reduce(nTotDissociationReactions, sumOp<label>());
            reduce(nTotIonisationReactions, sumOp<label>());
        }

        numberDensities_[0] = (molsReactants*cloud().nParticle())/volume;
        numberDensities_[1] = (molsReactants*cloud().nParticle())/volume;

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
            scalar reactionRateIonisation = 0.0;
            
            reactionRateDissociation =
            (
                nTotDissociationReactions
                * cloud_.nParticle()
            )/(counterIndex*deltaT*numberDensities_[0]* numberDensities_[1]*volume);
               
            Info<< "Dissociation type I reaction " 
                <<  reactantMolA << " + " << reactantMolB
                <<  " --> " 
                << productMolA << " + " << productMolB << " + " << reactantMolB 
                << ", reaction rate = " << reactionRateDissociation
                << endl;
                
            reactionRateIonisation =
            (
                nTotIonisationReactions
                * cloud_.nParticle()
            )/(counterIndex*deltaT*numberDensities_[0]* numberDensities_[1]*volume);
               
            Info<< "Ionisation type I reaction " 
                <<  reactantMolA << " + " << reactantMolB
                <<  " --> " 
                << productMolC << " + " << productMolD << " + " << reactantMolB 
                << ", reaction rate = " << reactionRateIonisation
                << endl;
        }
    }
    else
    {
        scalar nTotDissociationReactions = nDissociationReactions_;  
        scalar nTotIonisationReactions = nIonisationReactions_;  
        
        if(Pstream::parRun())
        {
            reduce(nTotDissociationReactions, sumOp<label>());
            reduce(nTotIonisationReactions, sumOp<label>());
        }
       
       if(nTotDissociationReactions > VSMALL)
       {
            word reactantMolA = cloud_.typeIdList()[reactantIds_[0]];
            word reactantMolB = cloud_.typeIdList()[reactantIds_[1]];

            word productMolA = cloud_.typeIdList()[productIdsDiss_[0]];
            word productMolB = cloud_.typeIdList()[productIdsDiss_[1]];
           
            Info<< "Dissociation type I reaction " 
                <<  reactantMolA << " + " << reactantMolB
                <<  " --> " 
                << productMolA << " + " << productMolB << " + " << reactantMolB 
                << " is active." << endl;
       }
       
       if(nTotIonisationReactions > VSMALL)
       {
            word reactantMolA = cloud_.typeIdList()[reactantIds_[0]];
            word reactantMolB = cloud_.typeIdList()[reactantIds_[1]];

            word productMolA = cloud_.typeIdList()[productIdsIon_[0]];
            word productMolB = cloud_.typeIdList()[productIdsIon_[1]];
           
            Info<< "Ionisation type I reaction " 
                <<  reactantMolA << " + " << reactantMolB
                <<  " --> " 
                << productMolA << " + " << productMolB << " + " << reactantMolB 
                << " is active." << endl;
       } 
    }

    nReactionsPerTimeStep_ = 0.0;

}


const bool& dissociationTypeISameSpecies::relax() const
{
    return relax_;
}

}
// End namespace Foam

// ************************************************************************* //
