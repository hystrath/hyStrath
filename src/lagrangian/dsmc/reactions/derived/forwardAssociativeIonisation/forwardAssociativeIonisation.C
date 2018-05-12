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

#include "forwardAssociativeIonisation.H"
#include "addToRunTimeSelectionTable.H"
#include "fvc.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(forwardAssociativeIonisation, 0);

addToRunTimeSelectionTable(dsmcReaction, forwardAssociativeIonisation, dictionary);



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
forwardAssociativeIonisation::forwardAssociativeIonisation
(
    Time& t,
    dsmcCloud& cloud,
    const dictionary& dict
)
:
    dsmcReaction(t, cloud, dict),
    propsDict_(dict.subDict(typeName + "Properties")),
    reactantIds_(),
    intermediateId_(),
    associativeIonisationProductIds_(),
    ionisationProductIds_(),
    reactionName_(propsDict_.lookup("reactionName")),
    heatOfReactionRecombination_(readScalar(propsDict_.lookup("heatOfReactionRecombination"))),
    heatOfReactionIntermediateIonisation_(readScalar(propsDict_.lookup("heatOfReactionIntermediateIonisation"))),
    heatOfReactionIonisation_(readScalar(propsDict_.lookup("heatOfReactionIonisation"))),
    nTotalIonisationReactions_(0),
    nTotalAssociativeIonisationReactions_(0),
    nIonisationReactionsPerTimeStep_(0),
    nAssociativeIonisationReactionsPerTimeStep_(0),
    relax_(true),
    allowSplitting_(false),
    writeRatesToTerminal_(false),
    nReactions_(0),
    volume_(0.0),
    numberDensities_(2, 0.0)
{

}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

forwardAssociativeIonisation::~forwardAssociativeIonisation()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void forwardAssociativeIonisation::initialConfiguration()
{
    setProperties();
}

void forwardAssociativeIonisation::setProperties()
{
    // reading in reactants

    const List<word> reactantMolecules (propsDict_.lookup("reactantAtoms"));

    if(reactantMolecules.size() != 2)
    {
        FatalErrorIn("forwardAssociativeIonisation::setProperties()")
            << "There should be two reactants atoms, instead of " 
            << reactantMolecules.size() << nl 
            << exit(FatalError);
    }
    
//     if(reactantMolecules[0] != reactantMolecules[1])
//     {
//         FatalErrorIn("forwardAssociativeIonisation::setProperties()")
//             << "Both reactant species must be the same, they are currently " 
// 	    << reactantMolecules[0] << " and " << reactantMolecules[1] << nl 
//             << exit(FatalError);
//     }

    reactantIds_.setSize(reactantMolecules.size(), -1);

    allowSplitting_ = Switch(propsDict_.lookup("allowSplitting"));
    
    writeRatesToTerminal_ = Switch(propsDict_.lookup("writeRatesToTerminal"));

    forAll(reactantIds_, r)
    {
        reactantIds_[r] = findIndex(cloud_.typeIdList(), reactantMolecules[r]);

        // check that reactants belong to the typeIdList (constant/dsmcProperties)
        if(reactantIds_[r] == -1)
        {
            FatalErrorIn("forwardAssociativeIonisation::setProperties()")
                << "Cannot find type id: " << reactantMolecules[r] << nl 
                << exit(FatalError);
        }

        // check that reactants are 'ATOMS' (not 'MOLECULES') 

        const scalar& rDofReactant = cloud_.constProps(reactantIds_[r]).rotationalDegreesOfFreedom();
    
        if(rDofReactant > VSMALL)
        {
            FatalErrorIn("forwardAssociativeIonisation::setProperties()")
                << "Reactant must be an atom (not a molecule): " << reactantMolecules[r] 
                << nl 
                << exit(FatalError);
        }
    }

    // reading in associative ionisation products

    const List<word> associativeIonisationProductMolecules (propsDict_.lookup("productsOfAssociativeIonisation"));
    
    associativeIonisationProductIds_.setSize(associativeIonisationProductMolecules.size(), -1);
    
    forAll(associativeIonisationProductIds_, i)
    {
        associativeIonisationProductIds_[i] = findIndex(cloud_.typeIdList(), associativeIonisationProductMolecules[i]);
        
        // check that products belong to the typeIdList (constant/dsmcProperties)
        if(associativeIonisationProductIds_[i] == -1)
        {
            FatalErrorIn("forwardAssociativeIonisation::setProperties()")
                << "Cannot find type id: " << associativeIonisationProductMolecules[i] << nl 
                << exit(FatalError);
        }
    }


    // check that products are a 'MOLECULE' and an 'ELECTRON'

    scalar rDofProd1 = cloud_.constProps(associativeIonisationProductIds_[0]).rotationalDegreesOfFreedom();

    if(rDofProd1 < 1)
    {
        FatalErrorIn("forwardAssociativeIonisation::setProperties()")
            << "First product must be a molecule: " << associativeIonisationProductMolecules 
            << nl 
            << exit(FatalError);
    }
    
    const label& charge = cloud_.constProps(associativeIonisationProductIds_[1]).charge();

    if(charge != -1)
    {
        FatalErrorIn("forwardAssociativeIonisation::setProperties()")
            << "Second product must be an electron: " << associativeIonisationProductMolecules 
            << nl 
            << exit(FatalError);
    }
    
    // reading in ionisation of P products

    const List<word> ionisationProductMolecules (propsDict_.lookup("productsOfIonisation"));
    
    ionisationProductIds_.setSize(ionisationProductMolecules.size(), -1);
    
    forAll(ionisationProductIds_, i)
    {
        ionisationProductIds_[i] = findIndex(cloud_.typeIdList(), ionisationProductMolecules[i]);
        
        // check that products belong to the typeIdList (constant/dsmcProperties)
        if(ionisationProductIds_[i] == -1)
        {
            FatalErrorIn("forwardAssociativeIonisation::setProperties()")
                << "Cannot find type id: " << ionisationProductMolecules[i] << nl 
                << exit(FatalError);
        }
    }


    // check that products are an 'ATOM' and an 'ELECTRON'

    rDofProd1 = cloud_.constProps(ionisationProductIds_[0]).rotationalDegreesOfFreedom();

    if(rDofProd1 > 0)
    {
        FatalErrorIn("forwardAssociativeIonisation::setProperties()")
            << "First product must be a charged atom: " << ionisationProductMolecules 
            << nl 
            << exit(FatalError);
    }
    
    const label& charge2 = cloud_.constProps(ionisationProductIds_[1]).charge();

    if(charge2 != -1)
    {
        FatalErrorIn("forwardAssociativeIonisation::setProperties()")
            << "Second product must be an electron: " << ionisationProductMolecules 
            << nl 
            << exit(FatalError);
    }
    
    //reading in intermediate molecule
    
    const word intermediateMolecule (propsDict_.lookup("intermediateMolecule"));
    
    intermediateId_ = findIndex(cloud_.typeIdList(), intermediateMolecule);

    // check that reactants belong to the typeIdList (constant/dsmcProperties)
    if(intermediateId_ == -1)
    {
        FatalErrorIn("forwardAssociativeIonisation::setProperties()")
            << "Cannot find type id: " << intermediateMolecule << nl 
            << exit(FatalError);
    }

    // check that the intermediate is a 'MOLECULE'

    const scalar& rDof = cloud_.constProps(intermediateId_).rotationalDegreesOfFreedom();

    if(rDof < 1)
    {
        FatalErrorIn("forwardAssociativeIonisation::setProperties()")
            << "The intermediate specie must be a molecule (not an atom): " << intermediateMolecule 
            << nl 
            << exit(FatalError);
    }

}

bool forwardAssociativeIonisation::tryReactMolecules(const label& typeIdP, const label& typeIdQ)
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

void forwardAssociativeIonisation::reaction
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


void forwardAssociativeIonisation::reaction
(
    dsmcParcel& p,
    dsmcParcel& q
)
{
    label typeIdP = p.typeId();
    label typeIdQ = q.typeId();
    
    if(reactantIds_[0] == reactantIds_[1] && typeIdP == reactantIds_[0]) //same species
    {
        relax_ = true;
        
        scalar totalReactionProbability = 0.0;
        scalarList reactionProbabilities(2, 0.0);
        
        vector UP = p.U();
        vector UQ = q.U();
        scalar mP = cloud_.constProps(typeIdP).mass();
        scalar mQ = cloud_.constProps(typeIdQ).mass();
        scalar mR = mP*mQ/(mP + mQ);
        scalar cRsqr = magSqr(UP - UQ);
        
        scalar translationalEnergy = 0.5*mR*cRsqr;
        
        scalar omegaIntermediate = cloud_.constProps(intermediateId_).omega();
        scalar rotationalDofIntermediate = cloud_.constProps(intermediateId_).rotationalDegreesOfFreedom();
        scalar ChiBIntermediate = 2.5 - omegaIntermediate;
        scalar thetaVIntermediate = cloud_.constProps(intermediateId_).thetaV()[0];
        scalar thetaDIntermediate = cloud_.constProps(intermediateId_).thetaD()[0];
        scalar ZrefIntermediate = cloud_.constProps(intermediateId_).Zref()[0];
        scalar refTempZvIntermediate = cloud_.constProps(intermediateId_).TrefZv()[0];
        label ELevelIntermediate = 0;
        List<scalar> EElistIntermediate = cloud_.constProps(intermediateId_).electronicEnergyList();
        List<label> gListIntermediate = cloud_.constProps(intermediateId_).degeneracyList();
        label jMaxIntermediate = cloud_.constProps(intermediateId_).numberOfElectronicLevels();
        
        scalar EEleP = cloud_.constProps(typeIdP).electronicEnergyList()[p.ELevel()];
        scalar EEleQ = cloud_.constProps(typeIdQ).electronicEnergyList()[q.ELevel()];
        
        label jMaxQ = cloud_.constProps(q.typeId()).numberOfElectronicLevels();
        List<scalar> EElistQ = cloud_.constProps(q.typeId()).electronicEnergyList();
        List<label> gListQ = cloud_.constProps(q.typeId()).degeneracyList();
        
        scalar omegaPQ =
        0.5
        *(
            cloud_.constProps(typeIdP).omega()
            + cloud_.constProps(typeIdQ).omega()
        );
        
        bool ionisationReaction = false;
        bool associativeIonisation = false;
                
        //2 reactions possible
        // 1. Ionisation of P
        // 2. Forward associative ionisation
        
        scalar Ec = 0.0;
        
        scalar ionisationEnergy = cloud_.constProps(typeIdP).ionisationTemperature()*physicoChemical::k.value();
        
        // calculate if an ionisation of species P is possible
        Ec = translationalEnergy + EEleP;

        if((Ec - ionisationEnergy) > VSMALL)
        {
            //Ionisation can occur
//             totalReactionProbability += 1.0;
//             reactionProbabilities[0] = 1.0;
        }
        
        //collision energy is the translational energy of the two atoms, plus their electronic energies
        
        Ec = translationalEnergy + EEleP;
        scalar EcOrig = Ec;  
         
        label iMax = (Ec /(physicoChemical::k.value()*thetaVIntermediate));
        
        label vibLevelIntermediate = -1;

        if(iMax > SMALL)
        {            
            vibLevelIntermediate = cloud_.postCollisionVibrationalEnergyLevel
                    (
                        true,
                        0.0,
                        iMax,
                        thetaVIntermediate,
                        thetaDIntermediate,
                        refTempZvIntermediate,
                        omegaIntermediate,
                        ZrefIntermediate,
                        Ec
                     );
        }
        
        if(vibLevelIntermediate == 0)
        {
            //'Form' the intermediate molecule and test it for ionisation
            const scalar& heatOfReactionRecombinationJoules = heatOfReactionRecombination_*physicoChemical::k.value();
            
            Ec = EcOrig + heatOfReactionRecombinationJoules;
            
            label postCollisionELevel = cloud_.postCollisionElectronicEnergyLevel
            (
                Ec,
                jMaxIntermediate,
                omegaIntermediate,
                EElistIntermediate,
                gListIntermediate
            );

            ELevelIntermediate = postCollisionELevel;
                
            // relative translational energy after electronic exchange
            Ec -=  EElistIntermediate[ELevelIntermediate];
            
            iMax = (Ec /(physicoChemical::k.value()*thetaVIntermediate));
            
            if(iMax > SMALL)
            {      
                label postCollisionVibLevel = cloud_.postCollisionVibrationalEnergyLevel
                    (
                        true,
                        0.0,
                        iMax,
                        thetaVIntermediate,
                        thetaDIntermediate,
                        refTempZvIntermediate,
                        omegaIntermediate,
                        ZrefIntermediate,
                        Ec
                     );

                Ec -= postCollisionVibLevel*thetaVIntermediate*physicoChemical::k.value();
            }
            
            scalar ERot = 0.0;
            
            scalar energyRatio = cloud_.postCollisionRotationalEnergy(rotationalDofIntermediate,ChiBIntermediate);

            ERot = energyRatio*Ec;
        
            Ec -= ERot;
    
            //redistribution finished, test it for ionisation
            
            scalar EcPPIon= 0.0;
            scalar ionisationEnergy = cloud_.constProps(intermediateId_).ionisationTemperature()*physicoChemical::k.value();
            
            // calculate if an ionisation of species P is possible
            EcPPIon = Ec + EElistIntermediate[ELevelIntermediate];
            
            if((EcPPIon - ionisationEnergy) > VSMALL)
            {
                //Associative ionisation can occur
                totalReactionProbability += 1.0;
                reactionProbabilities[1] = 1.0;
            }
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
                            ionisationReaction = true;
                            break;
                        }
                        if(i == 1)
                        {
                            //Associative ionisation reaction is to occur
                            associativeIonisation = true;
                            break;
                        }
                    }
                }
            }
        }
        
        if(ionisationReaction)
        {
            nTotalIonisationReactions_++;
            nIonisationReactionsPerTimeStep_++; 
            
            if(allowSplitting_)
            {
                relax_ = false;
                
                const scalar& heatOfReactionIonisationJoules = heatOfReactionIonisation_*physicoChemical::k.value();
                
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

                const label& typeId1 = ionisationProductIds_[0];
                const label& typeId2 = ionisationProductIds_[1];
                
                //Mass of Product one and two
                scalar mP1 = cloud_.constProps(typeId1).mass();
                scalar mP2 = cloud_.constProps(typeId2).mass();
                
                scalar mRatoms = mP1*mP2/(mP1 + mP2);
                
                translationalEnergy = 0.0;
                
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

                // Atom P will ionisie.
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
        
        if(associativeIonisation)
        {
            //IONISATION CAN OCCUR
            nTotalAssociativeIonisationReactions_++;
            nAssociativeIonisationReactionsPerTimeStep_++;
            
            if(allowSplitting_)
            {
                const scalar& heatOfReactionIonisationJoules = heatOfReactionIntermediateIonisation_*physicoChemical::k.value();
                const scalar& heatOfReactionRecombinationJoules = heatOfReactionRecombination_*physicoChemical::k.value();
                
                translationalEnergy += (heatOfReactionRecombinationJoules + heatOfReactionIonisationJoules);
                
                translationalEnergy += (EEleP + EEleQ);
                
                scalar thetaVNewP = cloud_.constProps(associativeIonisationProductIds_[0]).thetaV()[0];
                scalar thetaDNewP = cloud_.constProps(associativeIonisationProductIds_[0]).thetaD()[0];
                scalar jMaxNewP = cloud_.constProps(associativeIonisationProductIds_[0]).numberOfElectronicLevels();
                scalar rotationalDofNewP = cloud_.constProps(associativeIonisationProductIds_[0]).rotationalDegreesOfFreedom();
                scalar ZrefNewP = cloud_.constProps(associativeIonisationProductIds_[0]).Zref()[0];
                scalar refTempZvNewP = cloud_.constProps(associativeIonisationProductIds_[0]).TrefZv()[0];
                scalarList EElistNewP = cloud_.constProps(associativeIonisationProductIds_[0]).electronicEnergyList();
                labelList gListNewP = cloud_.constProps(associativeIonisationProductIds_[0]).degeneracyList();
                
                scalar omegaNewPQ =
                    0.5
                    *(
                        cloud_.constProps(associativeIonisationProductIds_[0]).omega()
                        + cloud_.constProps(associativeIonisationProductIds_[1]).omega()
                    );
                
                scalar ChiB = 2.5 - omegaNewPQ;
                
                label ELevelNewP = cloud_.postCollisionElectronicEnergyLevel
                                (
                                    translationalEnergy,
                                    jMaxNewP,
                                    omegaNewPQ,
                                    EElistNewP,
                                    gListNewP
                                );
                                
                translationalEnergy -= EElistNewP[ELevelNewP];
                
                label vibLevelNewP = 0;
                
                scalar ERotNewP = 0.0;

                if(rotationalDofNewP > VSMALL)
                {                    
                    label iMax = (translationalEnergy / (physicoChemical::k.value()*thetaVNewP));
                    
                    vibLevelNewP = cloud_.postCollisionVibrationalEnergyLevel
                                    (
                                            true,
                                            0,
                                            iMax,
                                            thetaVNewP,
                                            thetaDNewP,
                                            refTempZvNewP,
                                            omegaNewPQ,
                                            ZrefNewP,
                                            translationalEnergy
                                        );
                                    
                    translationalEnergy -= vibLevelNewP*thetaVNewP*physicoChemical::k.value();
                    
                    ERotNewP = translationalEnergy*cloud_.postCollisionRotationalEnergy(rotationalDofNewP,ChiB);
                            
                    translationalEnergy -= ERotNewP;
                }
                
                // centre of mass velocity of molecules (pre-split)
                vector Ucm = (mP*UP + mQ*UQ)/(mP + mQ);
                
                mP = cloud_.constProps(associativeIonisationProductIds_[0]).mass();
                mQ = cloud_.constProps(associativeIonisationProductIds_[1]).mass();
                
                mR = mP*mQ/(mP + mQ);
                
                scalar relVel = sqrt((2.0*translationalEnergy)/mR);

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
        
                UP = Ucm + (postCollisionRelU*mQ/(mP + mQ));
                UQ = Ucm - (postCollisionRelU*mP/(mP + mQ));
                
                p.typeId() = associativeIonisationProductIds_[0];
                p.U() = UP;
                p.ERot() = ERotNewP;
                p.vibLevel().setSize(1, vibLevelNewP);
                p.ELevel() = ELevelNewP;
                
                q.typeId() = associativeIonisationProductIds_[1];
                q.U() = UQ;
                q.ERot() = 0.0;
                q.vibLevel().setSize(0,0);
                q.ELevel() = 0;
            }
        }
    }
}

void  forwardAssociativeIonisation::outputResults(const label& counterIndex)
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
        label nTotalAssociativeIonisationReactions = nTotalAssociativeIonisationReactions_;
        label nTotalIonisationReactions = nTotalIonisationReactions_;

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
            reduce(nTotalAssociativeIonisationReactions, sumOp<label>());
            reduce(nTotalIonisationReactions, sumOp<label>());

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

        word productMolA = cloud_.typeIdList()[associativeIonisationProductIds_[0]];
        word productMolB = cloud_.typeIdList()[associativeIonisationProductIds_[1]];
        
        word productMolC = cloud_.typeIdList()[ionisationProductIds_[0]];
        word productMolD = cloud_.typeIdList()[ionisationProductIds_[1]];

        if((numberDensities_[0] > 0.0) && (numberDensities_[1] > 0.0))
        {
            scalar reactionRateAssociativeIonisation = 0.0;
            scalar reactionRateIonisation = 0.0;

            
            reactionRateAssociativeIonisation =
            (
                nTotalAssociativeIonisationReactions
                * cloud_.nParticle()
            )/(counterIndex*deltaT*numberDensities_[0]* numberDensities_[1]*volume);
            
            Info<< "Associative ionisation reaction " 
                <<  reactantMolA << " + " << reactantMolB
                <<  " --> " 
                << productMolA << " + " << productMolB
                << ", reaction rate = " << reactionRateAssociativeIonisation
                << endl;
                
            reactionRateIonisation =
            (
                nTotalIonisationReactions
                * cloud_.nParticle()
            )/(counterIndex*deltaT*numberDensities_[0]* numberDensities_[1]*volume);
            
            Info<< "Ionisation reaction " 
                <<  reactantMolA << " + " << reactantMolB
                <<  " --> " 
                << productMolC << " + " << productMolD << " + " << reactantMolB
                << ", reaction rate = " << reactionRateIonisation
                << endl;
        }
    }
    else
    {
        label nTotalAssociativeIonisationReactions = nTotalAssociativeIonisationReactions_;
        label nTotalIonisationReactions = nTotalIonisationReactions_;
        
        label nAssociativeIonisationReactionsPerTimeStep = nAssociativeIonisationReactionsPerTimeStep_;
        label nIonisationReactionsPerTimeStep = nIonisationReactionsPerTimeStep_;
        
        if(Pstream::parRun())
        {
            reduce(nTotalAssociativeIonisationReactions, sumOp<label>());
            reduce(nTotalIonisationReactions, sumOp<label>());
            
            reduce(nAssociativeIonisationReactionsPerTimeStep, sumOp<label>());
            reduce(nIonisationReactionsPerTimeStep, sumOp<label>());
        }
        
        if(nTotalAssociativeIonisationReactions > VSMALL)
        {
                word reactantMolA = cloud_.typeIdList()[reactantIds_[0]];
                word reactantMolB = cloud_.typeIdList()[reactantIds_[1]];

                word productMolA = cloud_.typeIdList()[associativeIonisationProductIds_[0]];
                word productMolB = cloud_.typeIdList()[associativeIonisationProductIds_[1]];

            
                Info<< "Associative ionisation reaction " 
                    <<  reactantMolA << " + " << reactantMolB
                    <<  " --> " 
                    << productMolA << " + " << productMolB
                    << " is active, nReactions this time step = " << nAssociativeIonisationReactionsPerTimeStep << endl;
        } 
        
        if(nTotalIonisationReactions > VSMALL)
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

    nAssociativeIonisationReactionsPerTimeStep_ = 0.0;
    nIonisationReactionsPerTimeStep_ = 0.0;
}


const bool& forwardAssociativeIonisation::relax() const
{
    return relax_;
}

}
// End namespace Foam

// ************************************************************************* //
