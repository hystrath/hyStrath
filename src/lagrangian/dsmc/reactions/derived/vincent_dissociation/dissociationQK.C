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

#include "dissociationQK.H"
#include "addToRunTimeSelectionTable.H"
#include "fvc.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(dissociationQK, 0);

addToRunTimeSelectionTable(dsmcReaction, dissociationQK, dictionary);


//

void dissociationQK::dissociateParticleByPartner(dsmcParcel& p, dsmcParcel& q)
{
    //dissociationQK of P only
    nTotDissociationReactions_++;
    nDissociationReactionsPerTimeStep_++;
    
    const label typeIdP = p.typeId();
    const label typeIdQ = q.typeId();
    
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
    
    const scalar thetaVQ = cloud_.constProps(typeIdQ).thetaV()[0];
    
    const scalar thetaDQ = cloud_.constProps(typeIdQ).thetaD()[0];
    
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
    
    if (allowSplitting_)
    {
        relax_ = false;
        
        const scalar omegaPQ =
            0.5
            *(
                    cloud_.constProps(typeIdP).omega()
                + cloud_.constProps(typeIdQ).omega()
            );
       
        const scalar ChiB = 2.5 - omegaPQ;
        
        translationalEnergy = translationalEnergy + heatOfReactionDissociationJoules_ + EVibP;
    
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
        
        ERotQ = translationalEnergy*cloud_.postCollisionRotationalEnergy(rotationalDofQ, ChiB);
                
        translationalEnergy -= ERotQ;
        
        scalar relVelNonDissoMol = sqrt(2.0*translationalEnergy/mR);

        //- Center of mass velocity of all particles
        vector Ucm = (mP*UP + mQ*UQ)/(mP + mQ);    

        //- Variable Hard Sphere collision part
        scalar cosTheta = 2.0*cloud_.rndGen().sample01<scalar>() - 1.0;
    
        scalar sinTheta = sqrt(1.0 - cosTheta*cosTheta);
    
        scalar phi = twoPi*cloud_.rndGen().sample01<scalar>();
    
        vector postCollisionRelU = relVelNonDissoMol
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
        
        //- Mass of Product one and two
        const scalar mP1 = cloud_.constProps(typeId1).mass();
        const scalar mP2 = cloud_.constProps(typeId2).mass();
        
        const scalar mRatoms = mP1*mP2/(mP1 + mP2);
        
        translationalEnergy = ERotP + EEleP;
        
        const scalar cRatoms = sqrt(2.0*translationalEnergy/mRatoms);

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
        label cell = -1;
        label tetFace = -1;
        label tetPt = -1;

        mesh_.findCellFacePt
        (
            p.position(),
            cell,
            tetFace,
            tetPt
        );
        
        p.typeId() = typeId1;
        p.U() = uP1;
        p.vibLevel().setSize(0,0);
        p.ERot() = 0.0;
        p.ELevel() = 0;
        
        labelList vibLevel;
        
        // insert new product 2
        cloud_.addNewParcel
        (
            p.position(),
            uP2,
            p.RWF(),
            0.0,
            0,
            cell,
            tetFace,
            tetPt,
            typeId2,
            -1,
            p.classification(),
            vibLevel
        );
    }
}





// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
dissociationQK::dissociationQK
(
    Time& t,
    dsmcCloud& cloud,
    const dictionary& dict
)
:
    dsmcReaction(t, cloud, dict),
    propsDict_(dict.subDict(typeName + "Properties")),
    reactantIds_(),
    productIdsDissociation_(),
    reactionName_(propsDict_.lookup("reactionName")),
    nTotDissociationReactions_(0),
    nDissociationReactionsPerTimeStep_(0),
    heatOfReactionDissociationJoules_
    (
        readScalar(propsDict_.lookup("heatOfReactionDissociation"))
       *physicoChemical::k.value()
    ),
    relax_(true),
    allowSplitting_(false),
    writeRatesToTerminal_(false),
    volume_(0.0),
    numberDensities_(2, 0.0)
{

}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

dissociationQK::~dissociationQK()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void dissociationQK::initialConfiguration()
{
    setProperties();
}


void dissociationQK::setProperties()
{
    allowSplitting_ = propsDict_.lookupOrDefault<Switch>("allowSplitting", true);
    
    writeRatesToTerminal_ = propsDict_.lookupOrDefault<Switch>("writeRatesToTerminal", false);
    
    //- Reading in reactants
    const List<word> reactantMolecules (propsDict_.lookup("reactantMolecules"));

    if (reactantMolecules.size() != 2)
    {
        FatalErrorIn("dissociationQK::setProperties()")
        << "There should be two reactants, instead of " 
        << reactantMolecules.size() << nl 
        << exit(FatalError);
    }
    
    if (reactantMolecules[0] != reactantMolecules[1])
    {
        FatalErrorIn("dissociationQK::setProperties()")
        << "Both reactant species must be the same, they are currently " 
	      << reactantMolecules[0] << " and " << reactantMolecules[1] << nl 
        << exit(FatalError);
    }

    reactantIds_.setSize(reactantMolecules.size(), -1);

    forAll(reactantIds_, r)
    {
        reactantIds_[r] = findIndex(cloud_.typeIdList(), reactantMolecules[r]);

        //- Check that reactants belong to the typeIdList (constant/dsmcProperties)
        if (reactantIds_[r] == -1)
        {
            FatalErrorIn("dissociationQK::setProperties()")
            << "Cannot find type id: " << reactantMolecules[r] << nl 
            << exit(FatalError);
        }

        //- Check that reactants are 'MOLECULES' (not 'ATOMS') 
        const label& rDof = cloud_.constProps(reactantIds_[r]).rotationalDegreesOfFreedom();
    
        if (rDof < 1)
        {
            FatalErrorIn("dissociationQK::setProperties()")
            << "Reactant must be a molecule (not an atom): " << reactantMolecules[r] 
            << nl 
            << exit(FatalError);
        }
        
        const label& vDof = cloud_.constProps(reactantIds_[r]).vibrationalDegreesOfFreedom();
        
        if (vDof > 1)
        {
            FatalErrorIn("dissociationQK::setProperties()")
            << "Reactions are currently only implemented for monatomic and diatomic species"
            << " This is a polyatomic:" << reactantMolecules[r] 
            << nl 
            << exit(FatalError);
        }
    }
    
    //- Reading in dissociationQK products
    const List<word> productMoleculesDissociation (propsDict_.lookup("productsOfDissociatedMolecule"));

    if (productMoleculesDissociation.size() != 2)
    {
        FatalErrorIn("dissociationQK::setProperties()")
        << "There should be two products, instead of " 
        << productMoleculesDissociation.size() << nl 
        << exit(FatalError);
    }
    
    productIdsDissociation_.setSize(productMoleculesDissociation.size(), -1);

    forAll(productIdsDissociation_, r)
    {
        productIdsDissociation_[r] = findIndex(cloud_.typeIdList(), productMoleculesDissociation[r]);

        //- Check that products belong to the typeIdList (constant/dsmcProperties)
        if (productIdsDissociation_[r] == -1)
        {
            FatalErrorIn("dissociationQK::setProperties()")
            << "Cannot find type id: " << productMoleculesDissociation[r] << nl 
            << exit(FatalError);
        }

        //- Check that products are 'ATOMS' (not 'MOLECULES') 
        const scalar& rDof = cloud_.constProps(productIdsDissociation_[r]).rotationalDegreesOfFreedom();
    
        if (rDof > 1)
        {
            FatalErrorIn("dissociationQK::setProperties()")
            << "Reactant must be an atom (not a molecule): " << productMoleculesDissociation[r] 
            << nl 
            << exit(FatalError);
        }
    }
}


bool dissociationQK::tryReactMolecules(const label& typeIdP, const label& typeIdQ)
{
    const label reactantPId = findIndex(reactantIds_, typeIdP);
    const label reactantQId = findIndex(reactantIds_, typeIdQ);

    if (reactantPId == reactantQId)
    {
        if ((reactantPId != -1) && (reactantQId != -1))
        {
            return true;
        }
    }

    if ((reactantPId != -1) && (reactantQId != -1) && (reactantPId != reactantQId))
    {
        return true;
    }
        
    return false;
}


void dissociationQK::reaction
(
    dsmcParcel& p,
    dsmcParcel& q,
    const DynamicList<label>& candidateList,
    const List<DynamicList<label> >& candidateSubList,
    const label& candidateP,
    const List<label>& whichSubCell
)
{}


void dissociationQK::reaction
(
    dsmcParcel& p,
    dsmcParcel& q
)
{
    const label typeIdP = p.typeId();
    const label typeIdQ = q.typeId();
    
    if (typeIdP == typeIdQ && typeIdP == reactantIds_[0]) // same species and desired species to measure rate for
    { 
        relax_ = true;
        
        scalar totalReactionProbability = 0.0;
        scalarList reactionProbabilities(2, 0.0);
        
        const scalar EVibP = p.vibLevel()[0]*cloud_.constProps(typeIdP).thetaV()[0]*physicoChemical::k.value();
        const scalar EVibQ = q.vibLevel()[0]*cloud_.constProps(typeIdQ).thetaV()[0]*physicoChemical::k.value();

        const scalar mP = cloud_.constProps(typeIdP).mass();
        const scalar mQ = cloud_.constProps(typeIdQ).mass();
        const scalar mR = mP*mQ/(mP + mQ);
        
        const scalar cRsqr = magSqr(p.U() - q.U());
        scalar translationalEnergy = 0.5*mR*cRsqr;
        
        const scalar thetaVP = cloud_.constProps(typeIdP).thetaV()[0];
        const scalar thetaVQ = cloud_.constProps(typeIdQ).thetaV()[0];
        
        const label idP = cloud_.constProps(typeIdP).charDissQuantumLevel()[0];
        const label idQ = cloud_.constProps(typeIdQ).charDissQuantumLevel()[0];

        const List<label>& gListP = cloud_.constProps(typeIdP).degeneracyList();
        const List<scalar>& EElistP = cloud_.constProps(typeIdP).electronicEnergyList();
        const List<label>& gListQ = cloud_.constProps(typeIdQ).degeneracyList();
        const List<scalar>& EElistQ = cloud_.constProps(typeIdQ).electronicEnergyList();

        //Test for P reactions first
        bool dissocReactionP = false;
                
        //1 reactions possible
        // 1. Dissociation of P

        const scalar EcPP = translationalEnergy + EVibP;
        const label imaxP = EcPP/(physicoChemical::k.value()*thetaVP);
        
        if (imaxP-idP > 0)
        {
            //Dissociation can occur
            totalReactionProbability += 1.0;
            reactionProbabilities[0] = 1.0;
        }
        
        //Decide if a reaction is to occur
        if (totalReactionProbability > cloud_.rndGen().sample01<scalar>())
        {
            //A chemical reaction is to occur, choose which one
            scalarList normalisedProbabilities(reactionProbabilities.size(), 0.0);
            scalar cumulativeProbability = 0.0;
            
            normalisedProbabilities = reactionProbabilities/totalReactionProbability;
            
            forAll(normalisedProbabilities, i)
            {                
                //If current reaction can't occur, don't check for it
                if (normalisedProbabilities[i] > VSMALL)
                {
                    cumulativeProbability += normalisedProbabilities[i];
                    
                    if (cumulativeProbability > cloud_.rndGen().sample01<scalar>())
                    {
                        //Current reaction is to occur
                        
                        if (i == 0)
                        {
                            //- Dissociation is to occur
                            dissocReactionP = true;
                            break;
                        }
                    }
                }
            }
        }
        
        //- Then test for Q reactions
        bool dissocReactionQ = false;
        
        scalar totalReactionProbabilityQ = 0.0;
        scalarList reactionProbabilitiesQ(1, 0.0);
        
        //1 reactions possible
        // 1. Dissociation of Q
        
        scalar EcPQ = translationalEnergy + EVibQ;
        
        if (dissocReactionP)
        {
            EcPQ += heatOfReactionDissociationJoules_;
        }
        
        label imaxQ = EcPQ/(physicoChemical::k.value()*thetaVQ);
        
        if (imaxQ-idQ > 0)
        {
            //Dissociation can occur
            totalReactionProbabilityQ += 1.0;
            reactionProbabilitiesQ[0] = 1.0;
        }
        
        //- Decide if a reaction is to occur
        if (totalReactionProbabilityQ > cloud_.rndGen().sample01<scalar>())
        {
            //- A chemical reaction is to occur, choose which one
            scalarList normalisedProbabilities(reactionProbabilitiesQ.size(), 0.0);
            scalar cumulativeProbability = 0.0;
            
            normalisedProbabilities = reactionProbabilitiesQ/totalReactionProbabilityQ;
            
            forAll(normalisedProbabilities, i)
            {                
                //If current reaction can't occur, don't check for it
                if (normalisedProbabilities[i] > VSMALL)
                {
                    cumulativeProbability += normalisedProbabilities[i];
                    
                    if (cumulativeProbability > cloud_.rndGen().sample01<scalar>())
                    {
                        // Dissociation reaction is to occur
                        if (i == 0)
                        {
                            //Dissociation is to occur
                            dissocReactionQ = true;
                            break;
                        }
                    }
                }
            }
        }
        
        
        //dissociationQK of P only
        if (dissocReactionP && !dissocReactionQ)
        {
            dissociateParticleByPartner(p, q);
        }
        
        //dissociationQK of Q only
        if (dissocReactionQ && !dissocReactionP)
        {
            dissociateParticleByPartner(q, p);
        }
        
        //dissociationQK of both P and Q 
        if (dissocReactionP && dissocReactionQ)
        {
            dissociateParticleByPartner(p, q);
            dissociateParticleByPartner(q, p);
        }
    }
}

void  dissociationQK::outputResults(const label& counterIndex)
{
    if (writeRatesToTerminal_ == true)
    {
        //- measure density 
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

                if (findIndex(reactantIds_, p->typeId()) != -1)
                {
                    molsReactants++;
                }
            }

            volume_ += mesh_.cellVolumes()[c];
        }
        
        scalar volume = volume_;
        label nTotReactions = nTotReactions_;

        //- Parallel communication
        if (Pstream::parRun())
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
        
        if ((numberDensities_[0] > 0.0) && (numberDensities_[1] > 0.0))
        {
            scalar reactionRateDissociation =
                (
                    nTotDissociationReactions_
                    * cloud_.nParticle()
                )
               /(counterIndex*deltaT*numberDensities_[0]* numberDensities_[1]*volume);
               
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
        label nDissociationReactionsPerTimeStep = nDissociationReactionsPerTimeStep_;
        
        if (Pstream::parRun())
        {
            reduce(nTotDissociationReactions, sumOp<label>());
            reduce(nDissociationReactionsPerTimeStep, sumOp<label>());
        }
       
       if (nTotDissociationReactions > VSMALL)
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
    }

    nDissociationReactionsPerTimeStep_ = 0;
}


const bool& dissociationQK::relax() const
{
    return relax_;
}

}
// End namespace Foam

// ************************************************************************* //
