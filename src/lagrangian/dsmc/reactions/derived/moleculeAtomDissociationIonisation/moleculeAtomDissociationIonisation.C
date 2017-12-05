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

#include "moleculeAtomDissociationIonisation.H"
#include "addToRunTimeSelectionTable.H"
#include "fvc.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(moleculeAtomDissociationIonisation, 0);

addToRunTimeSelectionTable(dsmcReaction, moleculeAtomDissociationIonisation, dictionary);



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
moleculeAtomDissociationIonisation::moleculeAtomDissociationIonisation
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
    productIdsIon2_(),
    reactionName_(propsDict_.lookup("reactionName")),
    heatOfReactionDiss_(readScalar(propsDict_.lookup("heatOfReactionDissociation"))),
    heatOfReactionIon_(readScalar(propsDict_.lookup("heatOfReactionIonisationMolecule"))),
    heatOfReactionIon2_(readScalar(propsDict_.lookup("heatOfReactionIonisationAtom"))),
    nTotReactionsDiss_(0),
    nTotReactionsIon_(0),
    nTotReactionsIon2_(0),
    nDissociationReactionsPerTimeStep_(0),
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

moleculeAtomDissociationIonisation::~moleculeAtomDissociationIonisation()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void moleculeAtomDissociationIonisation::initialConfiguration()
{
    setProperties();
}

void moleculeAtomDissociationIonisation::setProperties()
{
    // reading in reactants

    const List<word> reactantMolecules (propsDict_.lookup("reactants"));

    if(reactantMolecules.size() != 2)
    {
        FatalErrorIn("moleculeAtomDissociationIonisation::setProperties()")
            << "There should be two or more reactants, instead of " 
            << reactantMolecules.size() << nl 
            << exit(FatalError);
    }
    
    if(reactantMolecules[0] == reactantMolecules[1])
    {
        FatalErrorIn("moleculeAtomDissociationIonisation::setProperties()")
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
            FatalErrorIn("moleculeAtomDissociationIonisation::setProperties()")
                << "Cannot find type id: " << reactantMolecules[r] << nl 
                << exit(FatalError);
        }
    }
    
    // check that reactant one is a 'MOLECULE' (not an 'ATOM') 

    const label& rDof1 = cloud_.constProps(reactantIds_[0]).rotationalDegreesOfFreedom();

    if(rDof1 < 1)
    {
        FatalErrorIn("moleculeAtomDissociationIonisation::setProperties()")
            << "First reactant must be a molecule (not an atom): " << reactantMolecules[0] 
            << nl 
            << exit(FatalError);
    }
    
    const label& vDof1 = cloud_.constProps(reactantIds_[0]).vibrationalDegreesOfFreedom();

    if(vDof1 > 1)
    {
        FatalErrorIn("moleculeAtomDissociationIonisation::setProperties()")
            << "Reactions are currently only implemented for monatomic and diatomic species"
            << " This is a polyatomic:" << reactantMolecules[0] 
            << nl 
            << exit(FatalError);
    }
    
    // check that reactant two is an 'ATOM' (not n 'MOLECULE') 

    const label& rDof2 = cloud_.constProps(reactantIds_[1]).rotationalDegreesOfFreedom();

    if(rDof2 != 0)
    {
        FatalErrorIn("moleculeAtomDissociationIonisation::setProperties()")
            << "Second reactant must be an atom (not a molecule): " << reactantMolecules[1] 
            << nl 
            << exit(FatalError);
    }
    
    // reading in products

    List<word> productMoleculesDiss (propsDict_.lookup("productsOfDissociatedMolecule"));

    if(productMoleculesDiss.size() != 2)
    {
        FatalErrorIn("moleculeAtomDissociationIonisation::setProperties()")
            << "Number of dissociation products is " << productMoleculesDiss.size() <<
            ", should be two."
            << exit(FatalError);
    }
    
    List<word> productMoleculesIon (propsDict_.lookup("productsOfIonisedMolecule"));

    if(productMoleculesIon.size() != 2)
    {
        FatalErrorIn("moleculeAtomDissociationIonisation::setProperties()")
            << "Number of ionisation products is " << productMoleculesIon.size() <<
            ", should be two."
            << exit(FatalError);
    }
    
    List<word> productMoleculesIon2 (propsDict_.lookup("productsOfIonisedAtom"));

    if(productMoleculesIon2.size() != 2)
    {
        FatalErrorIn("moleculeAtomDissociationIonisation::setProperties()")
            << "Number of atom ionisation products is " << productMoleculesIon2.size() <<
            ", should be two."
            << exit(FatalError);
    }
    
/************************************************************************************************/

    productIdsDiss_.setSize(productMoleculesDiss.size());

    forAll(productMoleculesDiss, r)
    {
        if(productIdsDiss_.size() != 2)
        {
            FatalErrorIn("moleculeAtomDissociationIonisation::setProperties()")
                << "There should be two products (for the dissociating molecule "
                << reactantMolecules[r] << "), instead of " 
                << productIdsDiss_.size() << nl 
                << exit(FatalError);
        }
    
        forAll(productIdsDiss_, r)
        {
            productIdsDiss_[r] = findIndex(cloud_.typeIdList(), productMoleculesDiss[r]);

            // check that reactants belong to the typeIdList (constant/dsmcProperties)
            if(productIdsDiss_[r] == -1)
            {
                FatalErrorIn("moleculeAtomDissociationIonisation::setProperties()")
                    << "Cannot find type id: " << productMoleculesDiss[r] << nl 
                    << exit(FatalError);
            }
        }
        
        // check that product one is an 'ATOM' (not a 'MOLECULE') 

        const scalar& rDof3 = cloud_.constProps(productIdsDiss_[0]).rotationalDegreesOfFreedom();

        if(rDof3 != 0)
        {
            FatalErrorIn("moleculeAtomDissociationIonisation::setProperties()")
                << "First dissociation product must be an atom (not a molecule): " << productMoleculesDiss[0] 
                << nl 
                << exit(FatalError);
        }

        // check that product two is an 'ATOM' (not a 'MOLECULE') 

        const scalar& rDof4 = cloud_.constProps(productIdsDiss_[1]).rotationalDegreesOfFreedom();

        if(rDof4 != 0)
        {
            FatalErrorIn("moleculeAtomDissociationIonisation::setProperties()")
                << "Second dissociation product must be an atom (not a molecule): " << productMoleculesDiss[1] 
                << nl 
                << exit(FatalError);
        }
    }
    
/************************************************************************************************/

    productIdsIon_.setSize(productMoleculesIon.size());

    forAll(productMoleculesIon, r)
    {
        if(productIdsIon_.size() != 2)
        {
            FatalErrorIn("moleculeAtomDissociationIonisation::setProperties()")
                << "There should be two products (for the ionising molecule "
                << reactantMolecules[r] << "), instead of " 
                << productIdsIon_.size() << nl 
                << exit(FatalError);
        }
    
        forAll(productIdsIon_, r)
        {
            productIdsIon_[r] = findIndex(cloud_.typeIdList(), productMoleculesIon[r]);

            // check that reactants belong to the typeIdList (constant/dsmcProperties)
            if(productIdsIon_[r] == -1)
            {
                FatalErrorIn("moleculeAtomDissociationIonisation::setProperties()")
                    << "Cannot find type id: " << productMoleculesIon[r] << nl 
                    << exit(FatalError);
            }
        }
        
        // check that product one is an 'MOLECULE' (not an 'ATOM') 

        const scalar& rDof5 = cloud_.constProps(productIdsIon_[0]).rotationalDegreesOfFreedom();

        if(rDof5 < 0)
        {
            FatalErrorIn("moleculeAtomDissociationIonisation::setProperties()")
                << "First ionisation product must be a charged molecule (not an atom/electron): " << productMoleculesIon[0] 
                << nl 
                << exit(FatalError);
        }

        // check that product two is a 'ELECTRON'

        const label& charge = cloud_.constProps(productIdsIon_[1]).charge();

        if(charge != -1)
        {
            FatalErrorIn("moleculeAtomDissociationIonisation::setProperties()")
                << "Second ionisation product must be an electron: " << productMoleculesIon[1] 
                << nl 
                << exit(FatalError);
        }
    }
/************************************************************************************************/

    productIdsIon2_.setSize(productMoleculesIon2.size());

    forAll(productMoleculesIon2, r)
    {
        if(productIdsIon2_.size() != 2)
        {
            FatalErrorIn("moleculeAtomDissociationIonisation::setProperties()")
                << "There should be two products (for the ionising molecule "
                << reactantMolecules[r] << "), instead of " 
                << productIdsIon2_.size() << nl 
                << exit(FatalError);
        }
    
        forAll(productIdsIon2_, r)
        {
            productIdsIon2_[r] = findIndex(cloud_.typeIdList(), productMoleculesIon2[r]);

            // check that reactants belong to the typeIdList (constant/dsmcProperties)
            if(productIdsIon2_[r] == -1)
            {
                FatalErrorIn("moleculeAtomDissociationIonisation::setProperties()")
                    << "Cannot find type id: " << productMoleculesIon2[r] << nl 
                    << exit(FatalError);
            }
        }
        
        // check that product one is an 'ATOM'

        const scalar& rDof6 = cloud_.constProps(productIdsIon2_[0]).rotationalDegreesOfFreedom();

        if(rDof6 > 0)
        {
            FatalErrorIn("moleculeAtomDissociationIonisation::setProperties()")
                << "First ionisation product must be a charged atom: " << productMoleculesIon[0] 
                << nl 
                << exit(FatalError);
        }

        // check that product two is a 'ELECTRON'

        const label& charge = cloud_.constProps(productIdsIon2_[1]).charge();

        if(charge != -1)
        {
            FatalErrorIn("moleculeAtomDissociationIonisation::setProperties()")
                << "Second ionisation product must be an electron: " << productMoleculesIon[1] 
                << nl 
                << exit(FatalError);
        }
    }    
    
}

bool moleculeAtomDissociationIonisation::tryReactMolecules(const label& typeIdP, const label& typeIdQ)
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

void moleculeAtomDissociationIonisation::reaction
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


void moleculeAtomDissociationIonisation::reaction
(
    dsmcParcel& p,
    dsmcParcel& q
)
{
    label typeIdP = p.typeId();
    label typeIdQ = q.typeId();

    if(typeIdP == reactantIds_[0] && typeIdQ == reactantIds_[1]) // This produces the correct equilibrium rate A2 + X.
    {
        scalar totalReactionProbability = 0.0;
        scalarList reactionProbabilities(3, 0.0);
        
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
   
        scalar thetaVP = cloud_.constProps(typeIdP).thetaV()[0];
        List<label> gListP = cloud_.constProps(typeIdP).degeneracyList();
        List<scalar> EElistP = cloud_.constProps(typeIdP).electronicEnergyList();
        
        label jMaxP = cloud_.constProps(typeIdP).numberOfElectronicLevels();
        label jMaxQ = cloud_.constProps(typeIdQ).numberOfElectronicLevels();
        List<label> gListQ = cloud_.constProps(typeIdQ).degeneracyList();
        List<scalar> EElistQ = cloud_.constProps(typeIdQ).electronicEnergyList();

        scalar heatOfReactionJoulesDiss = heatOfReactionDiss_*physicoChemical::k.value();
        scalar heatOfReactionJoulesIon = heatOfReactionIon_*physicoChemical::k.value();
        scalar heatOfReactionJoulesIon2 = heatOfReactionIon2_*physicoChemical::k.value();

        scalar omegaPQ =
        0.5
        *(
            cloud_.constProps(typeIdP).omega()
            + cloud_.constProps(typeIdQ).omega()
        );
        
        bool dissocReaction = false;
        bool ionisationReaction = false;
        bool atomIonisationReaction = false;
        
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
        
        ionisationEnergy = cloud_.constProps(typeIdQ).ionisationTemperature()*physicoChemical::k.value();
        
        // calculate if an ionisation of species Q is possible
        EcPPIon = translationalEnergy + EEleQ;
        
        if((EcPPIon - ionisationEnergy) > VSMALL)
        {
            //IONISATION CAN OCCUR
            totalReactionProbability += 1.0;
            reactionProbabilities[1] = 1.0;
        }

        scalar EcPPDiss = 0.0;
        label idP = cloud_.constProps(typeIdP).charDissQuantumLevel()[0];
        label imaxP = 0;
        
        // calculate if a dissociation of species P is possible
        EcPPDiss = translationalEnergy + EVibP;

        imaxP = EcPPDiss/(physicoChemical::k.value()*thetaVP);

        if(imaxP-idP > 0)
        {
           //Dissociation can occur
            totalReactionProbability += 1.0;
            reactionProbabilities[2] = 1.0;
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
                            //Atom ionisation is to occur
                            atomIonisationReaction = true;
                            break;
                        }
                        if(i == 2)
                        {
                            //Dissociation reaction is to occur
                            dissocReaction = true;
                            break;
                        }
                    }
                }
            }
        }
        
        if(dissocReaction)
        {
            nDissociationReactionsPerTimeStep_++;
            nTotReactionsDiss_++;

            if(allowSplitting_)
            {
                relax_ = false;
                                
                translationalEnergy = translationalEnergy + heatOfReactionJoulesDiss + EVibP;
  
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
                UQ = Ucm - postCollisionRelU*mP/(mP + mQ); // Q is the NON-DISSOCIATING molecule.

                const label& typeId1 = productIdsDiss_[0];
                const label& typeId2 = productIdsDiss_[1];

                //Mass of Product one and two
                scalar mP1 = cloud_.constProps(typeId1).mass();
                scalar mP2 = cloud_.constProps(typeId2).mass();

                scalar mRatoms = mP1*mP2/(mP1 + mP2);
                
                //center of mass velocity of all particles

                vector UcmAtoms = UP;
                
                translationalEnergy = ERotP + EEleP;
                
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
        if(ionisationReaction)
        {
            nIonisationReactionsPerTimeStep_++;
            nTotReactionsIon_++;

            if(allowSplitting_)
            {
                relax_ = false;
                
                label ELevelQ = 0;
                                
                translationalEnergy = translationalEnergy + heatOfReactionJoulesIon + EEleP;

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
                scalar mP1 = cloud_.constProps(typeId1).mass(); // ion
                scalar mP2 = cloud_.constProps(typeId2).mass(); // electron

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
                p.vibLevel().setSize(1,0);
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
        if(atomIonisationReaction)
        {
            nIonisationReactions2PerTimeStep_++;
            nTotReactionsIon2_++;

            if(allowSplitting_)
            {
                relax_ = false;
                
                label ELevelP = 0;
                                
                translationalEnergy = translationalEnergy + heatOfReactionJoulesIon2 + EEleQ;
                
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
    
                UP = Ucm + postCollisionRelU*mQ/(mP + mQ); // P is the non-ionising molecule
                UQ = Ucm - postCollisionRelU*mP/(mP + mQ);

                const label& typeId1 = productIdsIon2_[0];
                const label& typeId2 = productIdsIon2_[1];
                
                //Mass of Product one and two
                scalar mP1 = cloud_.constProps(typeId1).mass(); // ion
                scalar mP2 = cloud_.constProps(typeId2).mass(); // electron

                scalar mRatoms = mP1*mP2/(mP1 + mP2);
                
                //center of mass velocity of all particles

                vector UcmAtoms = UQ;
                
                scalar translationalEnergy2 = 0;

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

                // New molecule P velocity.
                p.U() = UP;
                p.ELevel() = ELevelP;

                // Molecule P will ionise
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
  
    if(typeIdP == reactantIds_[1] && typeIdQ == reactantIds_[0]) // This produces the correct equilibrium rate X + A2.
    {
        scalar totalReactionProbability = 0.0;
        scalarList reactionProbabilities(3, 0.0);
        
        relax_ = true;
        
        vector UP = p.U();
        vector UQ = q.U();
        scalar ERotQ = q.ERot();
        scalar EVibQ = q.vibLevel()[0]*cloud_.constProps(typeIdQ).thetaV()[0]*physicoChemical::k.value();
        scalar EEleQ = cloud_.constProps(typeIdQ).electronicEnergyList()[q.ELevel()];
        scalar EEleP = cloud_.constProps(typeIdP).electronicEnergyList()[p.ELevel()];

        scalar mP = cloud_.constProps(typeIdP).mass();
        scalar mQ = cloud_.constProps(typeIdQ).mass();

        scalar mR = mP*mQ/(mP + mQ);
        scalar cRsqr = magSqr(UP - UQ);
        scalar translationalEnergy = 0.5*mR*cRsqr;
        
        scalar heatOfReactionJoulesDiss = heatOfReactionDiss_*physicoChemical::k.value();
        scalar heatOfReactionJoulesIon = heatOfReactionIon_*physicoChemical::k.value();
        scalar heatOfReactionJoulesIon2 = heatOfReactionIon2_*physicoChemical::k.value();
        
        label jMaxP = cloud_.constProps(typeIdP).numberOfElectronicLevels();
        label jMaxQ = cloud_.constProps(typeIdQ).numberOfElectronicLevels();
        List<label> gListP = cloud_.constProps(typeIdP).degeneracyList();
        List<scalar> EElistP = cloud_.constProps(typeIdP).electronicEnergyList();
        
        List<label> gListQ = cloud_.constProps(typeIdQ).degeneracyList();
        List<scalar> EElistQ = cloud_.constProps(typeIdQ).electronicEnergyList();

        scalar omegaPQ =
        0.5
        *(
            cloud_.constProps(typeIdP).omega()
            + cloud_.constProps(typeIdQ).omega()
        );
        
        bool dissocReaction = false;
        bool ionisationReaction = false;
        bool atomIonisationReaction = false;
        
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
        
        ionisationEnergy = cloud_.constProps(typeIdP).ionisationTemperature()*physicoChemical::k.value();
        
        // calculate if an ionisation of species P is possible
        EcPQIon = translationalEnergy + EEleP;
        
        if((EcPQIon - ionisationEnergy) > VSMALL)
        {
            //IONISATION CAN OCCUR
            totalReactionProbability += 1.0;
            reactionProbabilities[1] = 1.0;
        }
        
        scalar thetaVQ = cloud_.constProps(typeIdQ).thetaV()[0];

        scalar EcPQDiss = 0.0;
        label idQ = cloud_.constProps(typeIdQ).charDissQuantumLevel()[0];
        label imaxQ = 0;
        
        // calculate if a dissociation of species Q is possible
        EcPQDiss = translationalEnergy + EVibQ;

        imaxQ = EcPQDiss/(physicoChemical::k.value()*thetaVQ);

        if(imaxQ-idQ > 0)
        {
           //Dissociation can occur
            totalReactionProbability += 1.0;
            reactionProbabilities[2] = 1.0;
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
                            //Atom ionisation reaction is to occur
                            atomIonisationReaction = true;
                            break;
                        }
                        if(i == 2)
                        {
                            //Dissociation reaction is to occur
                            dissocReaction = true;
                            break;
                        }
                    }
                }
            }
        }

        if(dissocReaction)
        {
            nDissociationReactionsPerTimeStep_++;
            nTotReactionsDiss_++;
            
            if(allowSplitting_)
            {
                relax_ = false;
                
                translationalEnergy = translationalEnergy + heatOfReactionJoulesDiss + EVibQ;     
                
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
    
                UP = Ucm + (postCollisionRelU*mQ/(mP + mQ)); // UP is the single atom
                UQ = Ucm - (postCollisionRelU*mP/(mP + mQ)); // UQ is used as Ucm for atomic split.

                const label& typeId1 = productIdsDiss_[0];
                const label& typeId2 = productIdsDiss_[1];
                
                //Mass of Product one and two
                scalar mP1 = cloud_.constProps(typeId1).mass();
                scalar mP2 = cloud_.constProps(typeId2).mass();
                
                scalar mRatoms = mP1*mP2/(mP1 + mP2);
                
                translationalEnergy = ERotQ + EEleQ;
                
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


                vector uP1 = UQ + postCollisionRelU2*mP2/(mP1 + mP2);
                vector uP2 = UQ - postCollisionRelU2*mP1/(mP1 + mP2);

                // P remains NON-DISSOCIATED (it is an atom).
                p.U() = UP;
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
                q.U() = uP1;
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
                    uP2,
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
        
        if(ionisationReaction)
        {
            nIonisationReactionsPerTimeStep_++;
            nTotReactionsIon_++;

            if(allowSplitting_)
            {
                relax_ = false;
                
                label ELevelP = 0;
                                
                translationalEnergy = translationalEnergy + heatOfReactionJoulesIon + EEleQ;
                
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
                scalar mP1 = cloud_.constProps(typeId1).mass(); // ion
                scalar mP2 = cloud_.constProps(typeId2).mass(); // electron

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
                    0,
                    classificationP,
                    vibLevel
                );
            }
        }
        
        if(atomIonisationReaction)
        {
            nIonisationReactions2PerTimeStep_++;
            nTotReactionsIon2_++;

            if(allowSplitting_)
            {
                relax_ = false;
                
                label ELevelQ = 0;
                                
                translationalEnergy = translationalEnergy + heatOfReactionJoulesIon2 + EEleP;
                
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
                UQ = Ucm - postCollisionRelU*mP/(mP + mQ); 

                const label& typeId1 = productIdsIon2_[0];
                const label& typeId2 = productIdsIon2_[1];
                
                //Mass of Product one and two
                scalar mP1 = cloud_.constProps(typeId1).mass(); // ion
                scalar mP2 = cloud_.constProps(typeId2).mass(); // electron

                scalar mRatoms = mP1*mP2/(mP1 + mP2);
                
                //center of mass velocity of all particles

                vector UcmAtoms = UP;
                
                scalar translationalEnergy2 = 0;

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

                // New atom Q velocity.
                q.U() = UQ;
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
                    0,
                    classificationP,
                    vibLevel
                );
            }
        }
    }
}

void  moleculeAtomDissociationIonisation::outputResults(const label& counterIndex)
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
        label nTotReactionsDiss = nTotReactionsDiss_;
        label nTotReactionsIon = nTotReactionsIon_;
        label nTotReactionsIon2 = nTotReactionsIon2_;

        //- Parallel communication
        if(Pstream::parRun())
        {
            reduce(volume, sumOp<scalar>());
            reduce(mols[0], sumOp<label>());
            reduce(mols[1], sumOp<label>());
            reduce(nTotReactionsDiss, sumOp<label>());
            reduce(nTotReactionsIon, sumOp<label>());
            reduce(nTotReactionsIon2, sumOp<label>());
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
        
        word productMolE = cloud_.typeIdList()[productIdsIon2_[0]];
        word productMolF = cloud_.typeIdList()[productIdsIon2_[1]];

        if((numberDensities_[0] > 0.0) && (numberDensities_[1] > 0.0))
        { 

            scalar reactionRateDiss = 0.0;
            scalar reactionRateIon = 0.0;
            scalar reactionRateIon2 = 0.0;

            reactionRateDiss =
            (
                nTotReactionsDiss
                * cloud_.nParticle()
            )/(counterIndex*deltaT*numberDensities_[0]* numberDensities_[1]*volume);

            Info<< "Dissociation type II reaction "
                <<  reactantMolA << " + " << reactantMolB 
                <<  " --> "
                << productMolA << " + " << productMolB << " + " << reactantMolB 
                << ", reaction rate = " << reactionRateDiss
                << endl;
                
            reactionRateIon =
            (
                nTotReactionsIon
                * cloud_.nParticle()
            )/(counterIndex*deltaT*numberDensities_[0]* numberDensities_[1]*volume);

            Info<< "Ionisation reaction "
                <<  reactantMolA << " + " << reactantMolB 
                <<  " --> "
                << productMolC << " + " << productMolD << " + " << reactantMolB 
                << ", reaction rate = " << reactionRateIon
                << endl;
                
                reactionRateIon2 =
            (
                nTotReactionsIon2
                * cloud_.nParticle()
            )/(counterIndex*deltaT*numberDensities_[0]* numberDensities_[1]*volume);

            Info<< "Ionisation reaction "
                <<  reactantMolA << " + " << reactantMolB 
                <<  " --> "
                << reactantMolA << " + " << productMolE << " + " << productMolF 
                << ", reaction rate = " << reactionRateIon2
                << endl;
        }
    }
    else
    {
        label nTotReactionsDiss = nTotReactionsDiss_; 
        label nTotReactionsIon = nTotReactionsIon_; 
        label nTotReactionsIon2 = nTotReactionsIon2_; 
        
        label nDissociationReactionsPerTimeStep = nDissociationReactionsPerTimeStep_;
        label nIonisationReactionsPerTimeStep = nIonisationReactionsPerTimeStep_;
        label nIonisationReactions2PerTimeStep = nIonisationReactions2PerTimeStep_;
        
        if(Pstream::parRun())
        {
            reduce(nTotReactionsDiss, sumOp<label>());
            reduce(nTotReactionsIon, sumOp<label>());
            reduce(nTotReactionsIon2, sumOp<label>());
            
            reduce(nDissociationReactionsPerTimeStep, sumOp<label>());
            reduce(nIonisationReactionsPerTimeStep, sumOp<label>());
            reduce(nIonisationReactions2PerTimeStep, sumOp<label>());
        }
    
        if(nTotReactionsDiss > VSMALL)
        {
                word reactantMolA = cloud_.typeIdList()[reactantIds_[0]];
                word reactantMolB = cloud_.typeIdList()[reactantIds_[1]];

                word productMolA = cloud_.typeIdList()[productIdsDiss_[0]];
                word productMolB = cloud_.typeIdList()[productIdsDiss_[1]];
            
                Info<< "Dissociation type II reaction "
                <<  reactantMolA << " + " << reactantMolB 
                <<  " --> "
                << productMolA << " + " << productMolB << " + " << reactantMolB  
                    << " is active, nReactions this time step = " << nDissociationReactionsPerTimeStep << endl;
        }  
        
        if(nTotReactionsIon > VSMALL)
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
        
        if(nTotReactionsIon2 > VSMALL)
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

    nDissociationReactionsPerTimeStep_ = 0.0;
    nIonisationReactionsPerTimeStep_ = 0.0;
    nIonisationReactions2PerTimeStep_ = 0.0;

}


const bool& moleculeAtomDissociationIonisation::relax() const
{
    return relax_;
}

}
// End namespace Foam

// ************************************************************************* //
