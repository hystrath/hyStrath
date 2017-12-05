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

#include "atomAtomIonisationDissimilarSpecies.H"
#include "addToRunTimeSelectionTable.H"
#include "fvc.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(atomAtomIonisationDissimilarSpecies, 0);

addToRunTimeSelectionTable(dsmcReaction, atomAtomIonisationDissimilarSpecies, dictionary);



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
atomAtomIonisationDissimilarSpecies::atomAtomIonisationDissimilarSpecies
(
    Time& t,
    dsmcCloud& cloud,
    const dictionary& dict
)
:
    dsmcReaction(t, cloud, dict),
    propsDict_(dict.subDict(typeName + "Properties")),
    reactantIds_(),
    productIdsIon_(),
    reactionName_(propsDict_.lookup("reactionName")),
    chargedAtom_(false),
    heatOfReactionIon_(),
    heatOfReactionIon2_(),
    nTotIonisationReactions_(0),
    nTotIonisationReactions2_(0),
    nIonisationReactionsPerTimeStep_(0),
    nIonisationReactions2PerTimeStep_(0),
    relax_(true),
    allowSplitting_(false),
    writeRatesToTerminal_(false),
    volume_(0.0),
    numberDensities_(2, 0.0)
{

}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

atomAtomIonisationDissimilarSpecies::~atomAtomIonisationDissimilarSpecies()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void atomAtomIonisationDissimilarSpecies::initialConfiguration()
{
    setProperties();
}

void atomAtomIonisationDissimilarSpecies::setProperties()
{
    // reading in reactants

    const List<word> reactantMolecules (propsDict_.lookup("reactants"));

    if(reactantMolecules.size() != 2)
    {
        FatalErrorIn("atomAtomIonisationDissimilarSpecies::setProperties()")
            << "There should be two or more reactants, instead of " 
            << reactantMolecules.size() << nl 
            << exit(FatalError);
    }
    
    if(reactantMolecules[0] == reactantMolecules[1])
    {
        FatalErrorIn("atomAtomIonisationDissimilarSpecies::setProperties()")
            << "Reactant molecules cannot be same species." << nl
            << exit(FatalError);
    }

    reactantIds_.setSize(reactantMolecules.size(), -1);

    allowSplitting_ = Switch(propsDict_.lookup("allowSplitting"));
    
    chargedAtom_ = Switch(propsDict_.lookup("chargedAtom"));
    
    writeRatesToTerminal_ = Switch(propsDict_.lookup("writeRatesToTerminal"));

    forAll(reactantIds_, r)
    {
        reactantIds_[r] = findIndex(cloud_.typeIdList(), reactantMolecules[r]);

        // check that reactants belong to the typeIdList (constant/dsmcProperties)
        if(reactantIds_[r] == -1)
        {
            FatalErrorIn("atomAtomIonisationDissimilarSpecies::setProperties()")
                << "Cannot find type id: " << reactantMolecules[r] << nl 
                << exit(FatalError);
        }
    }
    
    // check that reactant one is an 'ATOM' 

    const scalar& rDof1 = cloud_.constProps(reactantIds_[0]).rotationalDegreesOfFreedom();

    if(rDof1 > VSMALL)
    {
        FatalErrorIn("atomAtomIonisationDissimilarSpecies::setProperties()")
            << "First reactant must be an atom (not a molecule or an electron): " << reactantMolecules[0] 
            << nl 
            << exit(FatalError);
    }
    
    // check that reactant two is an 'ATOM'

    const label& rDof2 = cloud_.constProps(reactantIds_[1]).rotationalDegreesOfFreedom();

    if(rDof2 > VSMALL)
    {
        FatalErrorIn("atomAtomIonisationDissimilarSpecies::setProperties()")
            << "Second reactant must be an atom (not a molecule or an electron): " << reactantMolecules[1] 
            << nl 
            << exit(FatalError);
    }
    
    const label& vDof1 = cloud_.constProps(reactantIds_[0]).vibrationalDegreesOfFreedom();

    if(vDof1 > VSMALL)
    {
         FatalErrorIn("atomIonIonisation::setProperties()")
            << "Reactions are currently only implemented for monatomic and diatomic species"
            << " This is a polyatomic:" << reactantMolecules[0] 
            << nl 
            << exit(FatalError);
    }
    
    // check that reactant two is an 'ATOM'

    const label& vDof2 = cloud_.constProps(reactantIds_[1]).vibrationalDegreesOfFreedom();

    if(vDof2 > VSMALL)
    {
         FatalErrorIn("atomIonIonisation::setProperties()")
            << "Reactions are currently only implemented for monatomic and diatomic species"
            << " This is a polyatomic:" << reactantMolecules[1] 
            << nl 
            << exit(FatalError);
    }
    
     const label& charge1 = cloud_.constProps(reactantIds_[0]).charge();

    if(charge1 == -1)
    {
        FatalErrorIn("atomAtomIonisationDissimilarSpecies::setProperties()")
            << "First reactant must be an atom (not a molecule or an electron): " << reactantMolecules[0] 
            << nl 
            << exit(FatalError);
    }
    
    // check that reactant two is an 'ATOM'

    const label& charge2 = cloud_.constProps(reactantIds_[1]).charge();

    if(charge2 == -1)
    {
        FatalErrorIn("atomAtomIonisationDissimilarSpecies::setProperties()")
            << "Second reactant must be an atom (not a molecule or an electron): " << reactantMolecules[1] 
            << nl 
            << exit(FatalError);
    }
    
    if(chargedAtom_)
    {
        // reading in ionisation products

        List<word> productMoleculesIonisation (propsDict_.lookup("productsOfIonisedAtom"));

        if(productMoleculesIonisation.size() != 2)
        {
            FatalErrorIn("atomAtomIonisationDissimilarSpecies::setProperties()")
                << "Number of ionisation products is " << productMoleculesIonisation.size() <<
                ", should be two."
                << exit(FatalError);
        }
        

        productIdsIon_.setSize(productMoleculesIonisation.size());

        forAll(productMoleculesIonisation, r)
        {
            if(productIdsIon_.size() != 2)
            {
                FatalErrorIn("atomAtomIonisationDissimilarSpecies::setProperties()")
                    << "There should be two products (for the ionising molecule "
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
                    FatalErrorIn("atomAtomIonisationDissimilarSpecies::setProperties()")
                        << "Cannot find type id: " << productMoleculesIonisation[r] << nl 
                        << exit(FatalError);
                }
            }
            
            // check that product one is a 'ATOM', not an 'MOLECULE'

            const scalar& rDof3 = cloud_.constProps(productIdsIon_[0]).rotationalDegreesOfFreedom();

            if(rDof3 > 1)
            {
                FatalErrorIn("atomAtomIonisationDissimilarSpecies::setProperties()")
                    << "First product must be an atom (not an atom): " << productMoleculesIonisation[0] 
                    << nl 
                    << exit(FatalError);
            }

            // check that product two is an 'ELECTRON'

            const label& charge = cloud_.constProps(productIdsIon_[1]).charge();

            if(charge != -1)
            {
                FatalErrorIn("atomAtomIonisationDissimilarSpecies::setProperties()")
                    << "Second product must be an electron: " << productMoleculesIonisation[1] 
                    << nl 
                    << exit(FatalError);
            }
        }
    }
    else
    {
        // reading in ionisation products

        List<word> productMoleculesIonisation (propsDict_.lookup("productsOfIonisedAtomP"));

        if(productMoleculesIonisation.size() != 2)
        {
            FatalErrorIn("atomAtomIonisationDissimilarSpecies::setProperties()")
                << "Number of ionisation products is " << productMoleculesIonisation.size() <<
                ", should be two."
                << exit(FatalError);
        }
        

        productIdsIon_.setSize(productMoleculesIonisation.size());

        forAll(productMoleculesIonisation, r)
        {
            if(productIdsIon_.size() != 2)
            {
                FatalErrorIn("atomAtomIonisationDissimilarSpecies::setProperties()")
                    << "There should be two products (for the ionising molecule "
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
                    FatalErrorIn("atomAtomIonisationDissimilarSpecies::setProperties()")
                        << "Cannot find type id: " << productMoleculesIonisation[r] << nl 
                        << exit(FatalError);
                }
            }
            
            // check that product one is a 'ATOM', not an 'MOLECULE'

            const scalar& rDof3 = cloud_.constProps(productIdsIon_[0]).rotationalDegreesOfFreedom();

            if(rDof3 > 1)
            {
                FatalErrorIn("atomAtomIonisationDissimilarSpecies::setProperties()")
                    << "First product must be an atom (not an atom): " << productMoleculesIonisation[0] 
                    << nl 
                    << exit(FatalError);
            }

            // check that product two is an 'ELECTRON'

            const label& charge = cloud_.constProps(productIdsIon_[1]).charge();

            if(charge != -1)
            {
                FatalErrorIn("atomAtomIonisationDissimilarSpecies::setProperties()")
                    << "Second product must be an electron: " << productMoleculesIonisation[1] 
                    << nl 
                    << exit(FatalError);
            }
        }
        
        List<word> productMoleculesIonisation2 (propsDict_.lookup("productsOfIonisedAtomQ"));

        if(productMoleculesIonisation2.size() != 2)
        {
            FatalErrorIn("atomAtomIonisationDissimilarSpecies::setProperties()")
                << "Number of ionisation products is " << productMoleculesIonisation2.size() <<
                ", should be two."
                << exit(FatalError);
        }
        

        productIdsIon2_.setSize(productMoleculesIonisation2.size());

        forAll(productMoleculesIonisation2, r)
        {
            if(productIdsIon2_.size() != 2)
            {
                FatalErrorIn("atomAtomIonisationDissimilarSpecies::setProperties()")
                    << "There should be two products (for the ionising molecule "
                    << reactantMolecules[r] << "), instead of " 
                    << productIdsIon2_.size() << nl 
                    << exit(FatalError);
            }
        
            forAll(productIdsIon2_, r)
            {
                productIdsIon2_[r] = findIndex(cloud_.typeIdList(), productMoleculesIonisation2[r]);

                // check that reactants belong to the typeIdList (constant/dsmcProperties)
                if(productIdsIon2_[r] == -1)
                {
                    FatalErrorIn("atomAtomIonisationDissimilarSpecies::setProperties()")
                        << "Cannot find type id: " << productMoleculesIonisation2[r] << nl 
                        << exit(FatalError);
                }
            }
            
            // check that product one is a 'ATOM', not an 'MOLECULE'

            const scalar& rDof3 = cloud_.constProps(productIdsIon2_[0]).rotationalDegreesOfFreedom();

            if(rDof3 > 1)
            {
                FatalErrorIn("atomAtomIonisationDissimilarSpecies::setProperties()")
                    << "First product must be an atom (not an atom): " << productMoleculesIonisation2[0] 
                    << nl 
                    << exit(FatalError);
            }

            // check that product two is an 'ELECTRON'

            const label& charge = cloud_.constProps(productIdsIon2_[1]).charge();

            if(charge != -1)
            {
                FatalErrorIn("atomAtomIonisationDissimilarSpecies::setProperties()")
                    << "Second product must be an electron: " << productMoleculesIonisation2[1] 
                    << nl 
                    << exit(FatalError);
            }
        }
    }

    if(chargedAtom_)
    {
        heatOfReactionIon_ = readScalar(propsDict_.lookup("heatOfReactionIonisation"));
    }
    else
    {
        heatOfReactionIon_ = readScalar(propsDict_.lookup("heatOfReactionIonisationP"));
        heatOfReactionIon2_ = readScalar(propsDict_.lookup("heatOfReactionIonisationQ"));
    }
}

bool atomAtomIonisationDissimilarSpecies::tryReactMolecules(const label& typeIdP, const label& typeIdQ)
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

void atomAtomIonisationDissimilarSpecies::reaction
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


void atomAtomIonisationDissimilarSpecies::reaction
(
    dsmcParcel& p,
    dsmcParcel& q
)
{
    label typeIdP = p.typeId();
    label typeIdQ = q.typeId();

    if( typeIdP == reactantIds_[0] && typeIdQ == reactantIds_[1] && !chargedAtom_ )
    {
        relax_ = true;
        
        scalar totalReactionProbability = 0.0;
        scalarList reactionProbabilities(2, 0.0);
        
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
        
        label jMaxP = cloud_.constProps(typeIdP).numberOfElectronicLevels();
        List<label> gListP = cloud_.constProps(typeIdP).degeneracyList();
        List<scalar> EElistP = cloud_.constProps(typeIdP).electronicEnergyList();
        
        label jMaxQ = cloud_.constProps(typeIdQ).numberOfElectronicLevels();
        List<label> gListQ = cloud_.constProps(typeIdQ).degeneracyList();
        List<scalar> EElistQ = cloud_.constProps(typeIdQ).electronicEnergyList();

        scalar heatOfReactionJoulesIon = heatOfReactionIon_*physicoChemical::k.value();
        scalar heatOfReactionJoulesIon2 = heatOfReactionIon2_*physicoChemical::k.value();
        
        scalar omegaPQ =
            0.5
            *(
                    cloud_.constProps(typeIdP).omega()
                + cloud_.constProps(typeIdQ).omega()
            );
            
        bool ionisationReactionP = false;
        bool ionisationReactionQ = false;
        
        scalar EcPPIon= 0.0;
        scalar ionisationEnergy = cloud_.constProps(typeIdP).ionisationTemperature()*physicoChemical::k.value();
        
        // calculate if an ionisation of species P is possible
        EcPPIon = translationalEnergy + EEleP;
        
        if((EcPPIon - ionisationEnergy) > VSMALL)
        {
            totalReactionProbability += 1.0;
            reactionProbabilities[0] = 1.0;
        }
        
        ionisationEnergy = cloud_.constProps(typeIdQ).ionisationTemperature()*physicoChemical::k.value();
        
        // calculate if an ionisation of species Q is possible
        EcPPIon = translationalEnergy + EEleQ;
        
        if((EcPPIon - ionisationEnergy) > VSMALL)
        {
            totalReactionProbability += 1.0;
            reactionProbabilities[1] = 1.0;
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
                            ionisationReactionP = true;
                            break;
                        }
                        if(i == 1)
                        {
                            //Ionisation reaction is to occur
                            ionisationReactionQ = true;
                            break;
                        }
                    }
                }
            }
        }
        
        if(ionisationReactionP)
        {
            nTotIonisationReactions_++;
            nIonisationReactionsPerTimeStep_++;
            
            if(allowSplitting_)
            {
                relax_ = false;
                                
                translationalEnergy = translationalEnergy + heatOfReactionJoulesIon + EEleP;
                
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
    
                UP = Ucm + postCollisionRelU*mQ/(mP + mQ);
                UQ = Ucm - postCollisionRelU*mP/(mP + mQ); // Q is the NON-IONISING atom.

                const label& typeId1 = productIdsIon_[0];
                const label& typeId2 = productIdsIon_[1];
                
                //Mass of Product one and two
                scalar mP1 = cloud_.constProps(typeId1).mass(); //
                scalar mP2 = cloud_.constProps(typeId2).mass(); //
                scalar mRatoms = mP1*mP2/(mP1 + mP2);
                
                //center of mass velocity of all particles

                vector UcmAtoms = UP;
                
                scalar translationalEnergy2 = ERotP + EVibP;

                scalar cRatoms = sqrt(2.0*translationalEnergy2/mRatoms);

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

                vector uP1 = UcmAtoms + postCollisionRelU2*mP2/(mP1 + mP2);
                vector uP2 = UcmAtoms - postCollisionRelU2*mP1/(mP1 + mP2);

                // New electron velocity.
                q.U() = UQ;
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
            nTotIonisationReactions2_++;
            nIonisationReactions2PerTimeStep_++;
            
            if(allowSplitting_)
            {
                relax_ = false;
                                
                translationalEnergy = translationalEnergy + heatOfReactionJoulesIon2 + EEleQ;
                
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
    
                UP = Ucm + postCollisionRelU*mQ/(mP + mQ); // P is the NON-IONISING atom.
                UQ = Ucm - postCollisionRelU*mP/(mP + mQ); // 

                const label& typeId1 = productIdsIon2_[0];
                const label& typeId2 = productIdsIon2_[1];
                
                //Mass of Product one and two
                scalar mP1 = cloud_.constProps(typeId1).mass(); //
                scalar mP2 = cloud_.constProps(typeId2).mass(); //

                scalar mRatoms = mP1*mP2/(mP1 + mP2);
                
                //center of mass velocity of all particles

                vector UcmAtoms = UQ;
                
                scalar translationalEnergy2 = ERotQ + EVibQ;

                scalar cRatoms = sqrt(2.0*translationalEnergy2/mRatoms);

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

                vector uQ1 = UcmAtoms + postCollisionRelU2*mP2/(mP1 + mP2);
                vector uQ2 = UcmAtoms - postCollisionRelU2*mP1/(mP1 + mP2);

                // New atom velocity.
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
    }
    
    if(typeIdP == reactantIds_[1] && typeIdQ == reactantIds_[0] && !chargedAtom_ ) // This produces the correct equilibrium rate A2 + X.
    {   
        relax_ = true;
        
        scalar totalReactionProbability = 0.0;
        scalarList reactionProbabilities(2, 0.0);
        
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
        
        label jMaxP = cloud_.constProps(typeIdP).numberOfElectronicLevels();
        List<label> gListP = cloud_.constProps(typeIdP).degeneracyList();
        List<scalar> EElistP = cloud_.constProps(typeIdP).electronicEnergyList();
        
        label jMaxQ = cloud_.constProps(typeIdQ).numberOfElectronicLevels();
        List<label> gListQ = cloud_.constProps(typeIdQ).degeneracyList();
        List<scalar> EElistQ = cloud_.constProps(typeIdQ).electronicEnergyList();

        scalar heatOfReactionJoulesIon = heatOfReactionIon_*physicoChemical::k.value();
        scalar heatOfReactionJoulesIon2 = heatOfReactionIon2_*physicoChemical::k.value();
        
        scalar omegaPQ =
            0.5
            *(
                    cloud_.constProps(typeIdP).omega()
                + cloud_.constProps(typeIdQ).omega()
            );
        
        bool ionisationReactionP = false;
        bool ionisationReactionQ = false;
        
        scalar EcPPIon= 0.0;
        scalar ionisationEnergy = cloud_.constProps(typeIdP).ionisationTemperature()*physicoChemical::k.value();
        
        // calculate if an ionisation of species P is possible
        EcPPIon = translationalEnergy + EEleP;
        
        if((EcPPIon - ionisationEnergy) > VSMALL)
        {
            totalReactionProbability += 1.0;
            reactionProbabilities[0] = 1.0;
        }
        
        ionisationEnergy = cloud_.constProps(typeIdQ).ionisationTemperature()*physicoChemical::k.value();
        
        // calculate if an ionisation of species Q is possible
        EcPPIon = translationalEnergy + EEleQ;
        
        if((EcPPIon - ionisationEnergy) > VSMALL)
        {
            totalReactionProbability += 1.0;
            reactionProbabilities[1] = 1.0;
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
                            ionisationReactionQ = true;
                            break;
                        }
                        if(i == 1)
                        {
                            //Ionisation reaction is to occur
                            ionisationReactionP = true;
                            break;
                        }
                    }
                }
            }
        }
        
        if(ionisationReactionP)
        {
            nTotIonisationReactions_++;
            nIonisationReactionsPerTimeStep_++;
            
            if(allowSplitting_)
            {
                relax_ = false;
                                
                translationalEnergy = translationalEnergy + heatOfReactionJoulesIon + EEleP;
                
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
    
                UP = Ucm + postCollisionRelU*mQ/(mP + mQ);
                UQ = Ucm - postCollisionRelU*mP/(mP + mQ); // Q is the NON-IONISING atom.

                const label& typeId1 = productIdsIon2_[0];
                const label& typeId2 = productIdsIon2_[1];
                
                //Mass of Product one and two
                scalar mP1 = cloud_.constProps(typeId1).mass(); //
                scalar mP2 = cloud_.constProps(typeId2).mass(); //
                scalar mRatoms = mP1*mP2/(mP1 + mP2);
                
                //center of mass velocity of all particles

                vector UcmAtoms = UP;
                
                scalar translationalEnergy2 = ERotP + EVibP;

                scalar cRatoms = sqrt(2.0*translationalEnergy2/mRatoms);

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

                vector uP1 = UcmAtoms + postCollisionRelU2*mP2/(mP1 + mP2);
                vector uP2 = UcmAtoms - postCollisionRelU2*mP1/(mP1 + mP2);

                // New electron velocity.
                q.U() = UQ;
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
            nTotIonisationReactions2_++;
            nIonisationReactions2PerTimeStep_++;
            
            if(allowSplitting_)
            {
                relax_ = false;
                                
                translationalEnergy = translationalEnergy + heatOfReactionJoulesIon2 + EEleQ;
                
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
    
                UP = Ucm + postCollisionRelU*mQ/(mP + mQ);
                UQ = Ucm - postCollisionRelU*mP/(mP + mQ); // Q is the NON-IONISING atom.

                const label& typeId1 = productIdsIon_[0];
                const label& typeId2 = productIdsIon_[1];
                
                //Mass of Product one and two
                scalar mP1 = cloud_.constProps(typeId1).mass(); //
                scalar mP2 = cloud_.constProps(typeId2).mass(); //
                scalar mRatoms = mP1*mP2/(mP1 + mP2);
                
                //center of mass velocity of all particles

                vector UcmAtoms = UQ;
                
                scalar translationalEnergy2 = ERotQ + EVibQ;

                scalar cRatoms = sqrt(2.0*translationalEnergy2/mRatoms);

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

                vector uQ1 = UcmAtoms + postCollisionRelU2*mP2/(mP1 + mP2);
                vector uQ2 = UcmAtoms - postCollisionRelU2*mP1/(mP1 + mP2);

                // New atom velocity.
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
        
    }
 /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 
    if( typeIdP == reactantIds_[0] && typeIdQ == reactantIds_[1] && chargedAtom_ )
    {
        relax_ = true;
        
        vector UP = p.U();
        vector UQ = q.U();
        scalar ERotP = p.ERot();
        scalar EVibP = p.vibLevel()[0]*cloud_.constProps(typeIdP).thetaV()[0]*physicoChemical::k.value();
        scalar EEleP = cloud_.constProps(typeIdP).electronicEnergyList()[p.ELevel()];
        scalar EEleQ = cloud_.constProps(typeIdQ).electronicEnergyList()[q.ELevel()];

        scalar mP = cloud_.constProps(typeIdP).mass();
        scalar mQ = cloud_.constProps(typeIdQ).mass();

        scalar mR = mP*mQ/(mP + mQ);
        scalar cRsqr = magSqr(UP - UQ);
        scalar translationalEnergy = 0.5*mR*cRsqr;
        
        List<label> gListP = cloud_.constProps(typeIdP).degeneracyList();
        List<scalar> EElistP = cloud_.constProps(typeIdP).electronicEnergyList();
        
        label jMaxQ = cloud_.constProps(typeIdQ).numberOfElectronicLevels();
        List<label> gListQ = cloud_.constProps(typeIdQ).degeneracyList();
        List<scalar> EElistQ = cloud_.constProps(typeIdQ).electronicEnergyList();

        scalar heatOfReactionJoulesIon = heatOfReactionIon_*physicoChemical::k.value();
        
        scalar omegaPQ =
            0.5
            *(
                    cloud_.constProps(typeIdP).omega()
                + cloud_.constProps(typeIdQ).omega()
            );
        
        scalar EcPPIon= 0.0;
        scalar ionisationEnergy = cloud_.constProps(typeIdP).ionisationTemperature()*physicoChemical::k.value();
        
        // calculate if an ionisation of species P is possible
        EcPPIon = translationalEnergy + EEleP;
        
        if((EcPPIon - ionisationEnergy) > VSMALL)
        {
            //P can ionise
            nTotIonisationReactions_++;
            nIonisationReactionsPerTimeStep_++;
            
            if(allowSplitting_)
            {
                relax_ = false;
                                
                translationalEnergy = translationalEnergy + heatOfReactionJoulesIon + EEleP;
                
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
    
                UP = Ucm + postCollisionRelU*mQ/(mP + mQ);
                UQ = Ucm - postCollisionRelU*mP/(mP + mQ); // Q is the NON-IONISING atom.

                const label& typeId1 = productIdsIon_[0];
                const label& typeId2 = productIdsIon_[1];
                
                //Mass of Product one and two
                scalar mP1 = cloud_.constProps(typeId1).mass(); //
                scalar mP2 = cloud_.constProps(typeId2).mass(); //
                scalar mRatoms = mP1*mP2/(mP1 + mP2);
                
                //center of mass velocity of all particles

                vector UcmAtoms = UP;
                
                scalar translationalEnergy2 = ERotP + EVibP;

                scalar cRatoms = sqrt(2.0*translationalEnergy2/mRatoms);

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

                vector uP1 = UcmAtoms + postCollisionRelU2*mP2/(mP1 + mP2);
                vector uP2 = UcmAtoms - postCollisionRelU2*mP1/(mP1 + mP2);

                // New electron velocity.
                q.U() = UQ;
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
    }
    
    if(typeIdP == reactantIds_[1] && typeIdQ == reactantIds_[0] && chargedAtom_ ) // This produces the correct equilibrium rate A2 + X.
    {   
        relax_ = true;
        
        vector UP = p.U();
        vector UQ = q.U();
        scalar ERotQ = q.ERot();
        scalar EVibQ = q.vibLevel()[0]*cloud_.constProps(typeIdQ).thetaV()[0]*physicoChemical::k.value();
        scalar EEleP = cloud_.constProps(typeIdP).electronicEnergyList()[p.ELevel()];
        scalar EEleQ = cloud_.constProps(typeIdQ).electronicEnergyList()[q.ELevel()];

        scalar mP = cloud_.constProps(typeIdP).mass();
        scalar mQ = cloud_.constProps(typeIdQ).mass();

        scalar mR = mP*mQ/(mP + mQ);
        scalar cRsqr = magSqr(UP - UQ);
        scalar translationalEnergy = 0.5*mR*cRsqr;
        
        label jMaxP = cloud_.constProps(typeIdP).numberOfElectronicLevels();
        List<label> gListP = cloud_.constProps(typeIdP).degeneracyList();
        List<scalar> EElistP = cloud_.constProps(typeIdP).electronicEnergyList();

        List<label> gListQ = cloud_.constProps(typeIdQ).degeneracyList();
        List<scalar> EElistQ = cloud_.constProps(typeIdQ).electronicEnergyList();

        scalar heatOfReactionJoulesIon = heatOfReactionIon_*physicoChemical::k.value();
        
        scalar omegaPQ =
            0.5
            *(
                    cloud_.constProps(typeIdP).omega()
                + cloud_.constProps(typeIdQ).omega()
            );
        
        // calculate if an ionisation of species Q is possible
        scalar EcPPIon = translationalEnergy + EEleQ;
        
        scalar ionisationEnergy = cloud_.constProps(typeIdQ).ionisationTemperature()*physicoChemical::k.value();
        
        if((EcPPIon - ionisationEnergy) > VSMALL)
        {
            //Q can ionise
            
            nTotIonisationReactions_++;
            nIonisationReactionsPerTimeStep_++;
            
            if(allowSplitting_)
            {
                relax_ = false;
                                
                translationalEnergy = translationalEnergy + heatOfReactionJoulesIon + EEleQ;
                
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
    
                UP = Ucm + postCollisionRelU*mQ/(mP + mQ);
                UQ = Ucm - postCollisionRelU*mP/(mP + mQ); // Q is the NON-IONISING atom.

                const label& typeId1 = productIdsIon_[0];
                const label& typeId2 = productIdsIon_[1];
                
                //Mass of Product one and two
                scalar mP1 = cloud_.constProps(typeId1).mass(); //
                scalar mP2 = cloud_.constProps(typeId2).mass(); //
                scalar mRatoms = mP1*mP2/(mP1 + mP2);
                
                //center of mass velocity of all particles

                vector UcmAtoms = UQ;
                
                scalar translationalEnergy2 = ERotQ + EVibQ;

                scalar cRatoms = sqrt(2.0*translationalEnergy2/mRatoms);

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

                vector uQ1 = UcmAtoms + postCollisionRelU2*mP2/(mP1 + mP2);
                vector uQ2 = UcmAtoms - postCollisionRelU2*mP1/(mP1 + mP2);

                // New atom velocity.
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
    }
}

void  atomAtomIonisationDissimilarSpecies::outputResults(const label& counterIndex)
{    
    if(writeRatesToTerminal_ == true)
    {
        const List< DynamicList<dsmcParcel*> >& cellOccupancy
            = cloud_.cellOccupancy();

        List<label> mols;
	mols.append(0);
	mols.append(0);
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
        label nTotReactionsIonisation = nTotIonisationReactions_;
        label nTotReactionsIonisation2 = nTotIonisationReactions2_;

        //- Parallel communication
        if(Pstream::parRun())
        {
            reduce(volume, sumOp<scalar>());
            reduce(mols[0], sumOp<label>());
            reduce(mols[1], sumOp<label>());
            reduce(nTotReactionsIonisation, sumOp<label>());
            reduce(nTotReactionsIonisation2, sumOp<label>());
        }

        numberDensities_[0] = (mols[0]*cloud().nParticle())/volume;
        numberDensities_[1] = (mols[1]*cloud().nParticle())/volume;

        const scalar& deltaT = mesh_.time().deltaT().value();

        word reactantMolA = cloud_.typeIdList()[reactantIds_[0]];
        word reactantMolB = cloud_.typeIdList()[reactantIds_[1]];
        
        word productMolC = cloud_.typeIdList()[productIdsIon_[0]];
        word productMolD = cloud_.typeIdList()[productIdsIon_[1]];
        
        word productMolE;
        word productMolF;
        
        if(!chargedAtom_)
        {
            productMolE = cloud_.typeIdList()[productIdsIon2_[0]];
            productMolF = cloud_.typeIdList()[productIdsIon2_[1]];
        }

        if((numberDensities_[0] > 0.0) && (numberDensities_[1] > 0.0))
        {   
            scalar reactionRateIonisation = 0.0;
            scalar reactionRateIonisation2 = 0.0;

            reactionRateIonisation =
            (
            nTotReactionsIonisation
            * cloud_.nParticle()
            )/(counterIndex*deltaT*numberDensities_[0]* numberDensities_[1]*volume);

            Info<< "Ionisation reaction "
                <<  reactantMolA << " + " << reactantMolB 
                <<  " --> "
                << productMolC << " + " << productMolD << " + " << reactantMolB 
                << ", reaction rate = " << reactionRateIonisation
                << endl;
                
            if(!chargedAtom_)
            {
                reactionRateIonisation2 =
                (
                nTotReactionsIonisation2
                * cloud_.nParticle()
                )/(counterIndex*deltaT*numberDensities_[0]* numberDensities_[1]*volume);

                Info<< "Ionisation reaction "
                    <<  reactantMolA << " + " << reactantMolB 
                    <<  " --> "
                    << reactantMolA << " + " << productMolE << " + " << productMolF 
                    << ", reaction rate = " << reactionRateIonisation2
                    << endl;
            }
        }
    }
    else
    {
        label nTotReactionsIonisation = nTotIonisationReactions_;   
        label nTotReactionsIonisation2 = nTotIonisationReactions2_;  
        label nIonisationReactionsPerTimeStep = nIonisationReactionsPerTimeStep_;
        label nIonisationReactions2PerTimeStep = nIonisationReactions2PerTimeStep_;
        
        if(Pstream::parRun())
        {
            reduce(nTotReactionsIonisation, sumOp<label>());
            reduce(nTotReactionsIonisation2, sumOp<label>());
            reduce(nIonisationReactionsPerTimeStep, sumOp<label>());
            reduce(nIonisationReactions2PerTimeStep, sumOp<label>());
        }  
        
        if(nTotReactionsIonisation > VSMALL)
        {
                word reactantMolA = cloud_.typeIdList()[reactantIds_[0]];
                word reactantMolB = cloud_.typeIdList()[reactantIds_[1]];

                word productMolA = cloud_.typeIdList()[productIdsIon_[0]];
                word productMolB = cloud_.typeIdList()[productIdsIon_[1]];
            
                Info<< "Ionisation reaction "
                <<  reactantMolA << " + " << reactantMolB 
                <<  " --> "
                << productMolA << " + " << productMolB << " + " << reactantMolB  
                    << " is active, nReactions this time step = " << nIonisationReactionsPerTimeStep << endl;
        } 
        
        if(nTotReactionsIonisation2 > VSMALL)
        {
                word reactantMolA = cloud_.typeIdList()[reactantIds_[0]];
                word reactantMolB = cloud_.typeIdList()[reactantIds_[1]];

                word productMolA = cloud_.typeIdList()[productIdsIon2_[0]];
                word productMolB = cloud_.typeIdList()[productIdsIon2_[1]];
            
                Info<< "Ionisation reaction "
                <<  reactantMolA << " + " << reactantMolB 
                <<  " --> "
                << reactantMolA << " + " << productMolA << " + " << productMolB  
                    << " is active, nReactions this time step = " << nIonisationReactions2PerTimeStep << endl;
        } 
    }

    nIonisationReactionsPerTimeStep_ = 0.0;
    nIonisationReactions2PerTimeStep_ = 0.0;
}


const bool& atomAtomIonisationDissimilarSpecies::relax() const
{
    return relax_;
}

}
// End namespace Foam

// ************************************************************************* //
