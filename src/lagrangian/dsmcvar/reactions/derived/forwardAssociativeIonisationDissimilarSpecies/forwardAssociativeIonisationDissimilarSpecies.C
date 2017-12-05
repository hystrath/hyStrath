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

#include "forwardAssociativeIonisationDissimilarSpecies.H"
#include "addToRunTimeSelectionTable.H"
#include "fvc.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(forwardAssociativeIonisationDissimilarSpecies, 0);

addToRunTimeSelectionTable(dsmcReaction, forwardAssociativeIonisationDissimilarSpecies, dictionary);



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
forwardAssociativeIonisationDissimilarSpecies::forwardAssociativeIonisationDissimilarSpecies
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
    ionisationPProductIds_(),
    ionisationQProductIds_(),
    reactionName_(propsDict_.lookup("reactionName")),
    heatOfReactionRecombination_(readScalar(propsDict_.lookup("heatOfReactionRecombination"))),
    heatOfReactionIntermediateIonisation_(readScalar(propsDict_.lookup("heatOfReactionIntermediateIonisation"))),
    heatOfReactionIonisationP_(readScalar(propsDict_.lookup("heatOfReactionIonisationP"))),
    heatOfReactionIonisationQ_(readScalar(propsDict_.lookup("heatOfReactionIonisationQ"))),
    nTotalIonisationPReactions_(0),
    nTotalIonisationQReactions_(0),
    nTotalAssociativeIonisationReactions_(0),
    nIonisationPReactionsPerTimeStep_(0),
    nIonisationQReactionsPerTimeStep_(0),
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

forwardAssociativeIonisationDissimilarSpecies::~forwardAssociativeIonisationDissimilarSpecies()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void forwardAssociativeIonisationDissimilarSpecies::initialConfiguration()
{
    setProperties();
}

void forwardAssociativeIonisationDissimilarSpecies::setProperties()
{
    // reading in reactants

    const List<word> reactantMolecules (propsDict_.lookup("reactantAtoms"));

    if(reactantMolecules.size() != 2)
    {
        FatalErrorIn("forwardAssociativeIonisationDissimilarSpecies::setProperties()")
            << "There should be two reactants atoms, instead of " 
            << reactantMolecules.size() << nl 
            << exit(FatalError);
    }
    
//     if(reactantMolecules[0] != reactantMolecules[1])
//     {
//         FatalErrorIn("forwardAssociativeIonisationDissimilarSpecies::setProperties()")
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
            FatalErrorIn("forwardAssociativeIonisationDissimilarSpecies::setProperties()")
                << "Cannot find type id: " << reactantMolecules[r] << nl 
                << exit(FatalError);
        }

        // check that reactants are 'ATOMS' (not 'MOLECULES') 

        const scalar& rDofReactant = cloud_.constProps(reactantIds_[r]).rotationalDegreesOfFreedom();
    
        if(rDofReactant > VSMALL)
        {
            FatalErrorIn("forwardAssociativeIonisationDissimilarSpecies::setProperties()")
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
            FatalErrorIn("forwardAssociativeIonisationDissimilarSpecies::setProperties()")
                << "Cannot find type id: " << associativeIonisationProductMolecules[i] << nl 
                << exit(FatalError);
        }
    }


    // check that products are a 'MOLECULE' and an 'ELECTRON'

    scalar rDofProd1 = cloud_.constProps(associativeIonisationProductIds_[0]).rotationalDegreesOfFreedom();

    if(rDofProd1 < 1)
    {
        FatalErrorIn("forwardAssociativeIonisationDissimilarSpecies::setProperties()")
            << "First product must be a molecule: " << associativeIonisationProductMolecules 
            << nl 
            << exit(FatalError);
    }
    
    const label& charge = cloud_.constProps(associativeIonisationProductIds_[1]).charge();

    if(charge != -1)
    {
        FatalErrorIn("forwardAssociativeIonisationDissimilarSpecies::setProperties()")
            << "Second product must be an electron: " << associativeIonisationProductMolecules 
            << nl 
            << exit(FatalError);
    }
    
    // reading in ionisation of P products

    const List<word> ionisationPProductMolecules (propsDict_.lookup("productsOfIonisationP"));
    
    ionisationPProductIds_.setSize(ionisationPProductMolecules.size(), -1);
    
    forAll(ionisationPProductIds_, i)
    {
        ionisationPProductIds_[i] = findIndex(cloud_.typeIdList(), ionisationPProductMolecules[i]);
        
        // check that products belong to the typeIdList (constant/dsmcProperties)
        if(ionisationPProductIds_[i] == -1)
        {
            FatalErrorIn("forwardAssociativeIonisationDissimilarSpecies::setProperties()")
                << "Cannot find type id: " << ionisationPProductMolecules[i] << nl 
                << exit(FatalError);
        }
    }


    // check that products are an 'ATOM' and an 'ELECTRON'

    rDofProd1 = cloud_.constProps(ionisationPProductIds_[0]).rotationalDegreesOfFreedom();

    if(rDofProd1 > 0)
    {
        FatalErrorIn("forwardAssociativeIonisationDissimilarSpecies::setProperties()")
            << "First product must be a charged atom: " << ionisationPProductMolecules 
            << nl 
            << exit(FatalError);
    }
    
    const label& charge2 = cloud_.constProps(ionisationPProductIds_[1]).charge();

    if(charge2 != -1)
    {
        FatalErrorIn("forwardAssociativeIonisationDissimilarSpecies::setProperties()")
            << "Second product must be an electron: " << ionisationPProductMolecules 
            << nl 
            << exit(FatalError);
    }
    
    // reading in ionisation of Q products

    const List<word> ionisationQProductMolecules (propsDict_.lookup("productsOfIonisationQ"));
    
    ionisationQProductIds_.setSize(ionisationQProductMolecules.size(), -1);
    
    forAll(ionisationQProductIds_, i)
    {
        ionisationQProductIds_[i] = findIndex(cloud_.typeIdList(), ionisationQProductMolecules[i]);
        
        // check that products belong to the typeIdList (constant/dsmcProperties)
        if(ionisationQProductIds_[i] == -1)
        {
            FatalErrorIn("forwardAssociativeIonisationDissimilarSpecies::setProperties()")
                << "Cannot find type id: " << ionisationQProductMolecules[i] << nl 
                << exit(FatalError);
        }
    }


    // check that products are an 'ATOM' and an 'ELECTRON'

    rDofProd1 = cloud_.constProps(ionisationQProductIds_[0]).rotationalDegreesOfFreedom();

    if(rDofProd1 > 0)
    {
        FatalErrorIn("forwardAssociativeIonisationDissimilarSpecies::setProperties()")
            << "First product must be a charged atom: " << ionisationQProductMolecules 
            << nl 
            << exit(FatalError);
    }
    
    const label& charge3 = cloud_.constProps(ionisationQProductIds_[1]).charge();

    if(charge3 != -1)
    {
        FatalErrorIn("forwardAssociativeIonisationDissimilarSpecies::setProperties()")
            << "Second product must be an electron: " << ionisationQProductMolecules 
            << nl 
            << exit(FatalError);
    }
    
    //reading in intermediate molecule
    
    const word intermediateMolecule (propsDict_.lookup("intermediateMolecule"));
    
    intermediateId_ = findIndex(cloud_.typeIdList(), intermediateMolecule);

    // check that reactants belong to the typeIdList (constant/dsmcProperties)
    if(intermediateId_ == -1)
    {
        FatalErrorIn("forwardAssociativeIonisationDissimilarSpecies::setProperties()")
            << "Cannot find type id: " << intermediateMolecule << nl 
            << exit(FatalError);
    }

    // check that the intermediate is a 'MOLECULE'

    const scalar& rDof = cloud_.constProps(intermediateId_).rotationalDegreesOfFreedom();

    if(rDof < 1)
    {
        FatalErrorIn("forwardAssociativeIonisationDissimilarSpecies::setProperties()")
            << "The intermediate specie must be a molecule (not an atom): " << intermediateMolecule 
            << nl 
            << exit(FatalError);
    }

}

bool forwardAssociativeIonisationDissimilarSpecies::tryReactMolecules(const label& typeIdP, const label& typeIdQ)
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

void forwardAssociativeIonisationDissimilarSpecies::reaction
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


void forwardAssociativeIonisationDissimilarSpecies::reaction
(
    dsmcParcel& p,
    dsmcParcel& q
)
{
    label typeIdP = p.typeId();
    label typeIdQ = q.typeId();

    if(typeIdP == reactantIds_[0] && typeIdQ == reactantIds_[1])
    {
        relax_ = true;
    
        scalar totalReactionProbability = 0.0;
        scalarList reactionProbabilities(3, 0.0);
    
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
        
        label jMaxP = cloud_.constProps(p.typeId()).numberOfElectronicLevels();
        List<scalar> EElistP = cloud_.constProps(p.typeId()).electronicEnergyList();
        List<label> gListP= cloud_.constProps(p.typeId()).degeneracyList();
        
        label jMaxQ = cloud_.constProps(q.typeId()).numberOfElectronicLevels();
        List<scalar> EElistQ = cloud_.constProps(q.typeId()).electronicEnergyList();
        List<label> gListQ = cloud_.constProps(q.typeId()).degeneracyList();
        
        scalar omegaPQ =
        0.5
        *(
            cloud_.constProps(typeIdP).omega()
            + cloud_.constProps(typeIdQ).omega()
        );
        
        //collision energy is the translational energy of the two atoms, plus their electronic energies
        
        bool ionisationReactionP = false;
        bool ionisationReactionQ = false;
        bool associativeIonisation = false;
                
        //3 reactions possible
        // 1. Ionisation of P
        // 2. Ionisation of Q
        // 3. Forward associative ionisation
        
        scalar Ec = 0.0;
        
        scalar ionisationEnergy = cloud_.constProps(typeIdP).ionisationTemperature()*physicoChemical::k.value();
        
        // calculate if an ionisation of species P is possible
        Ec = translationalEnergy + EEleP;

        if((Ec - ionisationEnergy) > VSMALL)
        {
            //Ionisation can occur
            totalReactionProbability += 1.0;
            reactionProbabilities[0] = 1.0;
        }
        
        ionisationEnergy = cloud_.constProps(typeIdQ).ionisationTemperature()*physicoChemical::k.value();
        
        // calculate if an ionisation of species P is possible
        Ec = translationalEnergy + EEleQ;

        if((Ec - ionisationEnergy) > VSMALL)
        {
            //Ionisation can occur
            totalReactionProbability += 1.0;
            reactionProbabilities[1] = 1.0;
        }
        
        Ec = translationalEnergy + EEleP;
        scalar EcOrig = Ec;        
        
        label iMax = (Ec /(physicoChemical::k.value()*thetaVIntermediate));
        
        label vibLevelIntermediate = -1;

        if(iMax > SMALL)
        {
            
            vibLevelIntermediate = cloud_.postCollisionVibrationalEnergyLevel
                (
                    false,
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
                reactionProbabilities[2] = 1.0;
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
                            //Ionisation is to occur
                            ionisationReactionP = true;
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
                            //Associative ionisation is to occur
                            associativeIonisation = true;
                            break;
                        }
                    }
                }
            }
        }
            
        if(ionisationReactionP)
        {
            nTotalIonisationPReactions_++;
            nIonisationPReactionsPerTimeStep_++;
            
            if(allowSplitting_)
            {
                relax_ = false;
                
                const scalar& heatOfReactionIonisationJoules = heatOfReactionIonisationP_*physicoChemical::k.value();
                
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

                const label& typeId1 = ionisationPProductIds_[0];
                const label& typeId2 = ionisationPProductIds_[1];
                
                //Mass of Product one and two
                scalar mP1 = cloud_.constProps(typeId1).mass();
                scalar mP2 = cloud_.constProps(typeId2).mass();
                
                scalar mRatoms = mP1*mP2/(mP1 + mP2);
                
                translationalEnergy = 0.0;
                
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
                    0,
                    classificationP,
                    vibLevel
                );
            }
        }
        
        if(ionisationReactionQ)
        {
            nTotalIonisationQReactions_++;
            nIonisationQReactionsPerTimeStep_++;
            
            if(allowSplitting_)
            {
                relax_ = false;
                
                const scalar& heatOfReactionIonisationJoules = heatOfReactionIonisationQ_*physicoChemical::k.value();
                
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

                const label& typeId1 = ionisationQProductIds_[0];
                const label& typeId2 = ionisationQProductIds_[1];
                
                //Mass of Product one and two
                scalar mP1 = cloud_.constProps(typeId1).mass();
                scalar mP2 = cloud_.constProps(typeId2).mass();
                
                scalar mRatoms = mP1*mP2/(mP1 + mP2);
                
                translationalEnergy = 0.0;
                
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

                // P remains NON-IONISED.
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
                    0,
                    classificationQ,
                    vibLevel
                );
            }
        }
        
        if(associativeIonisation)
        {
            //Associative IONISATION CAN OCCUR
            nTotalAssociativeIonisationReactions_++;
            nAssociativeIonisationReactionsPerTimeStep_++;
            
            if(allowSplitting_)
            {
                const scalar& heatOfReactionIonisationJoules = heatOfReactionIntermediateIonisation_*physicoChemical::k.value();
                const scalar& heatOfReactionRecombinationJoules = heatOfReactionRecombination_*physicoChemical::k.value();
                
                translationalEnergy += heatOfReactionRecombinationJoules + heatOfReactionIonisationJoules;
                
                translationalEnergy += EEleP + EEleQ;
                
                // centre of mass velocity of molecules (pre-split)
                vector Ucm = (mP*UP + mQ*UQ)/(mP + mQ);

                mP = cloud_.constProps(associativeIonisationProductIds_[0]).mass();
                mQ = cloud_.constProps(associativeIonisationProductIds_[1]).mass();
                
                mR = mP*mQ/(mP + mQ);
                
                scalar relVel = sqrt((2.0*translationalEnergy)/mR);

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
        
                UP = Ucm + (postCollisionRelU*mQ/(mP + mQ));
                UQ = Ucm - (postCollisionRelU*mP/(mP + mQ));
                
                p.typeId() = associativeIonisationProductIds_[0];
                p.U() = UP;
                p.ERot() = 0.0;
                p.vibLevel().setSize(1,0);
                p.ELevel() = 0;
                
                q.typeId() = associativeIonisationProductIds_[1];
                q.U() = UQ;
                q.ERot() = 0.0;
                q.vibLevel().setSize(0,0);
                q.ELevel() = 0;
            }
        }
    }
    
    if(typeIdP == reactantIds_[1] && typeIdQ == reactantIds_[0])
    {
        relax_ = true;
    
        scalar totalReactionProbability = 0.0;
        scalarList reactionProbabilities(3, 0.0);
    
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
        
        label jMaxP = cloud_.constProps(p.typeId()).numberOfElectronicLevels();
        List<scalar> EElistP = cloud_.constProps(p.typeId()).electronicEnergyList();
        List<label> gListP= cloud_.constProps(p.typeId()).degeneracyList();
        
        label jMaxQ = cloud_.constProps(q.typeId()).numberOfElectronicLevels();
        List<scalar> EElistQ = cloud_.constProps(q.typeId()).electronicEnergyList();
        List<label> gListQ = cloud_.constProps(q.typeId()).degeneracyList();
        
        scalar omegaPQ =
        0.5
        *(
            cloud_.constProps(typeIdP).omega()
            + cloud_.constProps(typeIdQ).omega()
        );
        
        bool ionisationReactionP = false; //N
        bool ionisationReactionQ = false; //O
        bool associativeIonisation = false;
                
        //3 reactions possible
        // 1. Ionisation of P
        // 2. Ionisation of Q
        // 3. Forward associative ionisation
        
        //collision energy is the translational energy of the two atoms, plus their electronic energies
        
        scalar Ec = 0.0;
        
        scalar ionisationEnergy = cloud_.constProps(typeIdP).ionisationTemperature()*physicoChemical::k.value();
        
        // calculate if an ionisation of species Q is possible
        Ec = translationalEnergy + EEleP;

        if((Ec - ionisationEnergy) > VSMALL)
        {
            //Ionisation can occur
            totalReactionProbability += 1.0;
            reactionProbabilities[1] = 1.0;
        }
        
        ionisationEnergy = cloud_.constProps(typeIdQ).ionisationTemperature()*physicoChemical::k.value();
        
        // calculate if an ionisation of species P is possible
        Ec = translationalEnergy + EEleQ;

        if((Ec - ionisationEnergy) > VSMALL)
        {
            //Ionisation can occur
            totalReactionProbability += 1.0;
            reactionProbabilities[0] = 1.0;
        }
        
        Ec = translationalEnergy + EEleQ;
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
                reactionProbabilities[2] = 1.0;
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
                            //Ionisation is to occur
                            ionisationReactionP = true;
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
                            //Associative ionisation is to occur
                            associativeIonisation = true;
                            break;
                        }
                    }
                }
            }
        }
            
        if(ionisationReactionP) //Q ionises, is called P here for consistency with the reaction rates
        {
            nTotalIonisationPReactions_++;
            nIonisationPReactionsPerTimeStep_++;
            
            if(allowSplitting_)
            {
                relax_ = false;
                
                const scalar& heatOfReactionIonisationJoules = heatOfReactionIonisationP_*physicoChemical::k.value();
                
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

                const label& typeId1 = ionisationPProductIds_[0];
                const label& typeId2 = ionisationPProductIds_[1];
                
                //Mass of Product one and two
                scalar mP1 = cloud_.constProps(typeId1).mass();
                scalar mP2 = cloud_.constProps(typeId2).mass();
                
                scalar mRatoms = mP1*mP2/(mP1 + mP2);
                
                translationalEnergy = 0.0;
                
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

                // P remains NON-IONISED.
                p.U() = UP;
                p.ELevel() = ELevelP;

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
                    0,
                    classificationQ,
                    vibLevel
                );
            }
        }
        
        if(ionisationReactionQ) //P ionises, is called Q here for consistency with the reaction rates
        {
            nTotalIonisationQReactions_++;
            nIonisationQReactionsPerTimeStep_++;
            
            if(allowSplitting_)
            {
                relax_ = false;
                
                const scalar& heatOfReactionIonisationJoules = heatOfReactionIonisationQ_*physicoChemical::k.value();
                
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
                                
                translationalEnergy -= EElistP[ELevelQ];
                
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

                const label& typeId1 = ionisationQProductIds_[0];
                const label& typeId2 = ionisationQProductIds_[1];
                
                //Mass of Product one and two
                scalar mP1 = cloud_.constProps(typeId1).mass();
                scalar mP2 = cloud_.constProps(typeId2).mass();
                
                scalar mRatoms = mP1*mP2/(mP1 + mP2);
                
                translationalEnergy = 0.0;
                
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

                // Q remains NON-IONISED.
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
                    0,
                    classificationP,
                    vibLevel
                );
            }
        }
        
        if(associativeIonisation)
        {
            nTotalAssociativeIonisationReactions_++;
            nAssociativeIonisationReactionsPerTimeStep_++;
            
            if(allowSplitting_)
            {
                const scalar& heatOfReactionIonisationJoules = heatOfReactionIntermediateIonisation_*physicoChemical::k.value();
                const scalar& heatOfReactionRecombinationJoules = heatOfReactionRecombination_*physicoChemical::k.value();
                
                translationalEnergy += heatOfReactionRecombinationJoules + heatOfReactionIonisationJoules;
                
                translationalEnergy += EEleP + EEleQ;
                
                // centre of mass velocity of molecules (pre-split)
                vector Ucm = (mP*UP + mQ*UQ)/(mP + mQ);

                mP = cloud_.constProps(associativeIonisationProductIds_[0]).mass();
                mQ = cloud_.constProps(associativeIonisationProductIds_[1]).mass();
                
                mR = mP*mQ/(mP + mQ);
                
                scalar relVel = sqrt((2.0*translationalEnergy)/mR);

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
        
                UP = Ucm + (postCollisionRelU*mQ/(mP + mQ));
                UQ = Ucm - (postCollisionRelU*mP/(mP + mQ));
                
                p.typeId() = associativeIonisationProductIds_[0];
                p.U() = UP;
                p.ERot() = 0.0;
                p.vibLevel().setSize(1,0);
                p.ELevel() = 0;
                
                q.typeId() = associativeIonisationProductIds_[1];
                q.U() = UQ;
                q.ERot() = 0.0;
                q.vibLevel().setSize(0,0);
                q.ELevel() = 0;
            }
        }
    }
}

void  forwardAssociativeIonisationDissimilarSpecies::outputResults(const label& counterIndex)
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
	mols.append(0);
	mols.append(0);
        scalar volume = volume_;
        label nTotalAssociativeIonisationReactions = nTotalAssociativeIonisationReactions_;
        label nTotalIonisationPReactions = nTotalIonisationPReactions_;
        label nTotalIonisationQReactions = nTotalIonisationQReactions_;

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
            reduce(nTotalIonisationPReactions, sumOp<label>());
            reduce(nTotalIonisationQReactions, sumOp<label>());
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
        
        word productMolC = cloud_.typeIdList()[ionisationPProductIds_[0]];
        word productMolD = cloud_.typeIdList()[ionisationPProductIds_[1]];
        
        word productMolE = cloud_.typeIdList()[ionisationQProductIds_[0]];
        word productMolF = cloud_.typeIdList()[ionisationQProductIds_[1]];

        if((numberDensities_[0] > 0.0) && (numberDensities_[1] > 0.0))
        {
            scalar reactionRateAssociativeIonisation = 0.0;
            scalar reactionRateIonisationP = 0.0;
            scalar reactionRateIonisationQ = 0.0;
            
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
                
            reactionRateIonisationP =
            (
                nTotalIonisationPReactions
                * cloud_.nParticle()
            )/(counterIndex*deltaT*numberDensities_[0]* numberDensities_[1]*volume);
            
            Info<< "Ionisation reaction " 
                <<  reactantMolA << " + " << reactantMolB
                <<  " --> " 
                << productMolC << " + " << productMolD << " + " << reactantMolB
                << ", reaction rate = " << reactionRateIonisationP
                << endl;
                
            reactionRateIonisationQ =
            (
                nTotalIonisationQReactions
                * cloud_.nParticle()
            )/(counterIndex*deltaT*numberDensities_[0]* numberDensities_[1]*volume);
            
            Info<< "Ionisation reaction " 
                <<  reactantMolA << " + " << reactantMolB
                <<  " --> " 
                << reactantMolA << " + " << productMolE << " + " << productMolF
                << ", reaction rate = " << reactionRateIonisationQ
                << endl;
        }
    }
    else
    {
        label nTotalAssociativeIonisationReactions = nTotalAssociativeIonisationReactions_;
        label nTotalIonisationPReactions = nTotalIonisationPReactions_;
        label nTotalIonisationQReactions = nTotalIonisationQReactions_;
        
        label nAssociativeIonisationReactionsPerTimeStep = nAssociativeIonisationReactionsPerTimeStep_;
        label nIonisationPReactionsPerTimeStep = nIonisationPReactionsPerTimeStep_;
        label nIonisationQReactionsPerTimeStep = nIonisationQReactionsPerTimeStep_;
        
        if(Pstream::parRun())
        {
            reduce(nTotalAssociativeIonisationReactions, sumOp<label>());
            reduce(nTotalIonisationPReactions, sumOp<label>());
            reduce(nTotalIonisationQReactions, sumOp<label>());
            
            reduce(nAssociativeIonisationReactionsPerTimeStep, sumOp<label>());
            reduce(nIonisationPReactionsPerTimeStep, sumOp<label>());
            reduce(nIonisationQReactionsPerTimeStep, sumOp<label>());
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
        
        if(nTotalIonisationPReactions > VSMALL)
        {
                word reactantMolA = cloud_.typeIdList()[reactantIds_[0]];
                word reactantMolB = cloud_.typeIdList()[reactantIds_[1]];

                word productMolA = cloud_.typeIdList()[ionisationPProductIds_[0]];
                word productMolB = cloud_.typeIdList()[ionisationPProductIds_[1]];

            
                Info<< "Ionisation reaction " 
                    <<  reactantMolA << " + " << reactantMolB
                    <<  " --> " 
                    << productMolA << " + " << productMolB << " + " << reactantMolB
                    << " is active, nReactions this time step = " << nIonisationPReactionsPerTimeStep << endl;
        } 
        
        if(nTotalIonisationQReactions > VSMALL)
        {
                word reactantMolA = cloud_.typeIdList()[reactantIds_[0]];
                word reactantMolB = cloud_.typeIdList()[reactantIds_[1]];

                word productMolA = cloud_.typeIdList()[ionisationQProductIds_[0]];
                word productMolB = cloud_.typeIdList()[ionisationQProductIds_[1]];

            
                Info<< "Ionisation reaction " 
                    <<  reactantMolA << " + " << reactantMolB
                    <<  " --> " 
                    << reactantMolA << " + " << productMolA << " + " << productMolB
                    << " is active, nReactions this time step = " << nIonisationQReactionsPerTimeStep << endl;
        } 
    }

    nAssociativeIonisationReactionsPerTimeStep_ = 0.0;
    nIonisationPReactionsPerTimeStep_ = 0.0;
    nIonisationQReactionsPerTimeStep_ = 0.0;
}


const bool& forwardAssociativeIonisationDissimilarSpecies::relax() const
{
    return relax_;
}

}
// End namespace Foam

// ************************************************************************* //
