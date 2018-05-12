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

#include "reverseAssociativeIonisation.H"
#include "addToRunTimeSelectionTable.H"
#include "fvc.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(reverseAssociativeIonisation, 0);

addToRunTimeSelectionTable(dsmcReaction, reverseAssociativeIonisation, dictionary);



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
reverseAssociativeIonisation::reverseAssociativeIonisation
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
    productIds_(),
    reactionName_(propsDict_.lookup("reactionName")),
    heatOfReactionDissociation_(readScalar(propsDict_.lookup("heatOfReactionDissociation"))),
    heatOfReactionRecombination_(readScalar(propsDict_.lookup("heatOfReactionRecombination"))),
    relax_(true),
    allowSplitting_(false),
    writeRatesToTerminal_(false),
    nReactions_(0),
    volume_(0.0),
    numberDensities_(2, 0.0)
{

}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

reverseAssociativeIonisation::~reverseAssociativeIonisation()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void reverseAssociativeIonisation::initialConfiguration()
{
    setProperties();
}

void reverseAssociativeIonisation::setProperties()
{
    // reading in reactants

    const List<word> reactantMolecules (propsDict_.lookup("reactants"));

    if(reactantMolecules.size() != 2)
    {
        FatalErrorIn("reverseAssociativeIonisation::setProperties()")
            << "There should be two reactants, instead of " 
            << reactantMolecules.size() << nl 
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
            FatalErrorIn("reverseAssociativeIonisation::setProperties()")
                << "Cannot find type id: " << reactantMolecules[r] << nl 
                << exit(FatalError);
        }
    }
    
    // check that first reactant is a charged molecule

    const label& rDofReactant1 = cloud_.constProps(reactantIds_[0]).rotationalDegreesOfFreedom();
    const label& chargeReactant1 = cloud_.constProps(reactantIds_[0]).charge();
    
    if(rDofReactant1 < VSMALL )
    {
        FatalErrorIn("reverseAssociativeIonisation::setProperties()")
            << "First reactant must be an ionised molecule: " << reactantMolecules[0] 
            << nl 
            << exit(FatalError);
    }
    
    if(chargeReactant1 != 1 )
    {
        FatalErrorIn("reverseAssociativeIonisation::setProperties()")
            << "First reactant must be an ionised molecule: " << reactantMolecules[0] 
            << nl 
            << exit(FatalError);
    }
    
    const scalar& vDofReactant1 = cloud_.constProps(reactantIds_[0]).vibrationalDegreesOfFreedom();

    if(vDofReactant1 > 1)
    {
        FatalErrorIn("reverseAssociativeIonisation::setProperties()")
            << "Reactions are currently only implemented for monatomic and diatomic species"
            << " This is a polyatomic:" << reactantMolecules[0] 
            << nl 
            << exit(FatalError);
    }
    
    // check that second reactant is an electron

    const label& charge = cloud_.constProps(reactantIds_[1]).charge();

    if(charge != -1)
    {
        FatalErrorIn("reverseAssociativeIonisation::setProperties()")
            << "Second reactant must be an electron: " << reactantMolecules[0] 
            << nl 
            << exit(FatalError);
    }

    // reading in product

    const List<word> productMolecules (propsDict_.lookup("productsOfAssociativeIonisation"));
    
    productIds_.setSize(productMolecules.size(), -1);
    
    forAll(productIds_, i)
    {
        productIds_[i] = findIndex(cloud_.typeIdList(), productMolecules[i]);
        
        // check that products belong to the typeIdList (constant/dsmcProperties)
        if(productIds_[i] == -1)
        {
            FatalErrorIn("reverseAssociativeIonisation::setProperties()")
                << "Cannot find type id: " << productMolecules[i] << nl 
                << exit(FatalError);
        }
        
            // check that products are 'ATOMS' and not 'MOLECULES'

        const scalar& rDof = cloud_.constProps(productIds_[i]).rotationalDegreesOfFreedom();

        if(rDof > 1)
        {
            FatalErrorIn("reverseAssociativeIonisation::setProperties()")
                << "Product must be an atom: " << productMolecules[i] 
                << nl 
                << exit(FatalError);
        }
    }

    
    //reading in intermediate molecule
    
    const word intermediateMolecule (propsDict_.lookup("intermediateMolecule"));
    
    intermediateId_ = findIndex(cloud_.typeIdList(), intermediateMolecule);

    // check that reactants belong to the typeIdList (constant/dsmcProperties)
    if(intermediateId_ == -1)
    {
        FatalErrorIn("reverseAssociativeIonisation::setProperties()")
            << "Cannot find type id: " << intermediateMolecule << nl 
            << exit(FatalError);
    }

    // check that the intermediate is a 'MOLECULE'

    const scalar& rDof = cloud_.constProps(intermediateId_).rotationalDegreesOfFreedom();

    if(rDof < 1)
    {
        FatalErrorIn("reverseAssociativeIonisation::setProperties()")
            << "The intermediate specie must be a molecule (not an atom): " << intermediateMolecule 
            << nl 
            << exit(FatalError);
    }

}

bool reverseAssociativeIonisation::tryReactMolecules(const label& typeIdP, const label& typeIdQ)
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

void reverseAssociativeIonisation::reaction
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


void reverseAssociativeIonisation::reaction
(
    dsmcParcel& p,
    dsmcParcel& q
)
{
    label typeIdP = p.typeId();
    label typeIdQ = q.typeId();
    
    if(typeIdP == reactantIds_[0] && typeIdQ == reactantIds_[1])
    {
        //Perform a 'recombination' to test the intermediate molecule
        
        relax_ = true;
    
        vector UP = p.U();
        vector UQ = q.U();
        scalar mP = cloud_.constProps(typeIdP).mass();
        scalar mQ = cloud_.constProps(typeIdQ).mass();
        scalar mR = mP*mQ/(mP + mQ);
        scalar cRsqr = magSqr(UP - UQ);
        
        scalar translationalEnergy = 0.5*mR*cRsqr;
        
        scalar ERotP = p.ERot();
        scalar EVibP = p.vibLevel()[0]*cloud_.constProps(p.typeId()).thetaV()[0]*physicoChemical::k.value();
        
        scalar omegaIntermediate = cloud_.constProps(intermediateId_).omega();
        scalar rotationalDofIntermediate = cloud_.constProps(intermediateId_).rotationalDegreesOfFreedom();
        scalar ChiBIntermediate = 2.5 - omegaIntermediate;
        scalar thetaVIntermediate = cloud_.constProps(intermediateId_).thetaV()[0];
        scalar thetaDIntermediate = cloud_.constProps(intermediateId_).thetaD()[0];
        scalar ZrefIntermediate = cloud_.constProps(intermediateId_).Zref()[0];
        scalar refTempZvIntermediate = cloud_.constProps(intermediateId_).TrefZv()[0];
        label ELevelIntermediate = -1;
        List<scalar> EElistIntermediate = cloud_.constProps(intermediateId_).electronicEnergyList();
        List<label> gListIntermediate = cloud_.constProps(intermediateId_).degeneracyList();
        label jMaxIntermediate = cloud_.constProps(intermediateId_).numberOfElectronicLevels();
        
        scalar EEleP = cloud_.constProps(typeIdP).electronicEnergyList()[p.ELevel()];
        scalar EEleQ = cloud_.constProps(typeIdQ).electronicEnergyList()[q.ELevel()];
        
        //collision energy is the translational energy of the two atoms, plus their electronic energies
        
        scalar Ec = translationalEnergy + EEleP;
        scalar EcOrig = Ec;        
        
        label iMax = (Ec /(physicoChemical::k.value()*thetaVIntermediate));

        label postCollisionELevel = cloud_.postCollisionElectronicEnergyLevel
                            (
                                Ec,
                                jMaxIntermediate,
                                omegaIntermediate,
                                EElistIntermediate,
                                gListIntermediate
                            );

        ELevelIntermediate = postCollisionELevel;
        
        if(ELevelIntermediate == 0)
        {
            //'Form' the intermediate molecule and test it for ionisation
            const scalar& heatOfReactionRecombinationJoules = heatOfReactionRecombination_*physicoChemical::k.value();
            
            Ec = EcOrig + heatOfReactionRecombinationJoules;
            
            ELevelIntermediate = cloud_.postCollisionElectronicEnergyLevel
                            (
                                Ec,
                                jMaxIntermediate,
                                omegaIntermediate,
                                EElistIntermediate,
                                gListIntermediate
                            );

            iMax = (Ec /(physicoChemical::k.value()*thetaVIntermediate));
            
            label eVibLevel = -1;
            
            if(iMax > SMALL)
            {      
                eVibLevel = cloud_.postCollisionVibrationalEnergyLevel
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

                Ec -= eVibLevel*physicoChemical::k.value()*thetaVIntermediate;
            }
                
            // relative translational energy after electronic exchange
            Ec -=  EElistIntermediate[ELevelIntermediate];
            
            scalar ERot = 0.0;
            
            scalar energyRatio = cloud_.postCollisionRotationalEnergy(rotationalDofIntermediate,ChiBIntermediate);

            ERot = energyRatio*Ec;
        
            Ec -= ERot;

            //redistribution finished, test it for dissociation
            
            scalar EcDiss = 0.0;
            label iMax = 0;
            label id = cloud_.constProps(intermediateId_).charDissQuantumLevel()[0];
            
            // calculate if a dissociation is possible
            scalar EVib = eVibLevel*physicoChemical::k.value()*thetaVIntermediate;
            
            EcDiss = Ec + EVib;
            iMax = EcDiss/(physicoChemical::k.value()*thetaVIntermediate);
            
            if((iMax - id) > VSMALL)
            {
                //ASS. ION. CAN OCCUR
                nReactions_++;
                nReactionsPerTimeStep_++;
                
                if(allowSplitting_)
                {
                    const scalar& heatOfReactionDissociation = heatOfReactionDissociation_*physicoChemical::k.value();
                    
                    translationalEnergy += heatOfReactionDissociation + heatOfReactionRecombinationJoules;
                    
                    translationalEnergy += EEleP + EEleQ + EVibP + ERotP;
                    
                    // centre of mass velocity of molecules (pre-split)
                    vector Ucm = (mP*UP + mQ*UQ)/(mP + mQ);
                    
                    mP = cloud_.constProps(productIds_[0]).mass();
                    mQ = cloud_.constProps(productIds_[1]).mass();
                    
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
                    
                    p.typeId() = productIds_[1];
                    p.U() = UP;
                    p.ERot() = 0.0;
                    p.vibLevel().setSize(0,0);
                    p.ELevel() = 0;
                    
                    q.typeId() = productIds_[0];
                    q.U() = UQ;
                    q.ERot() = 0.0;
                    q.vibLevel().setSize(0,0);
                    q.ELevel() = 0;
                }
            }  
        }
    }
    
    if(typeIdP == reactantIds_[1] && typeIdQ == reactantIds_[0])
    {
        //Perform a 'recombination' to test the intermediate molecule
        
        relax_ = true;
    
        vector UP = p.U();
        vector UQ = q.U();
        scalar mP = cloud_.constProps(typeIdP).mass();
        scalar mQ = cloud_.constProps(typeIdQ).mass();
        scalar mR = mP*mQ/(mP + mQ);
        scalar cRsqr = magSqr(UP - UQ);
        
        scalar translationalEnergy = 0.5*mR*cRsqr;

        scalar ERotQ = q.ERot();
        scalar EVibQ = q.vibLevel()[0]*cloud_.constProps(q.typeId()).thetaV()[0]*physicoChemical::k.value();
        
        scalar omegaIntermediate = cloud_.constProps(intermediateId_).omega();
        scalar rotationalDofIntermediate = cloud_.constProps(intermediateId_).rotationalDegreesOfFreedom();
        scalar ChiBIntermediate = 2.5 - omegaIntermediate;
        scalar thetaVIntermediate = cloud_.constProps(intermediateId_).thetaV()[0];
        scalar thetaDIntermediate = cloud_.constProps(intermediateId_).thetaD()[0];
        scalar ZrefIntermediate = cloud_.constProps(intermediateId_).Zref()[0];
        scalar refTempZvIntermediate = cloud_.constProps(intermediateId_).TrefZv()[0];
        label ELevelIntermediate = -1;
        List<scalar> EElistIntermediate = cloud_.constProps(intermediateId_).electronicEnergyList();
        List<label> gListIntermediate = cloud_.constProps(intermediateId_).degeneracyList();
        label jMaxIntermediate = cloud_.constProps(intermediateId_).numberOfElectronicLevels();
        
        scalar EEleP = cloud_.constProps(typeIdP).electronicEnergyList()[p.ELevel()];
        scalar EEleQ = cloud_.constProps(typeIdQ).electronicEnergyList()[q.ELevel()];
        
        //collision energy is the translational energy of the two atoms, plus their electronic energies
        
        scalar Ec = translationalEnergy + EEleQ;
        scalar EcOrig = Ec;        
        
        label iMax = (Ec /(physicoChemical::k.value()*thetaVIntermediate));

        label postCollisionELevel = cloud_.postCollisionElectronicEnergyLevel
                            (
                                Ec,
                                jMaxIntermediate,
                                omegaIntermediate,
                                EElistIntermediate,
                                gListIntermediate
                            );

        ELevelIntermediate = postCollisionELevel;
        
        if(ELevelIntermediate == 0)
        {
            //'Form' the intermediate molecule and test it for ionisation
            const scalar& heatOfReactionRecombinationJoules = heatOfReactionRecombination_*physicoChemical::k.value();
            
            Ec = EcOrig + heatOfReactionRecombinationJoules;
            
            ELevelIntermediate = cloud_.postCollisionElectronicEnergyLevel
                            (
                                Ec,
                                jMaxIntermediate,
                                omegaIntermediate,
                                EElistIntermediate,
                                gListIntermediate
                            );
            
            iMax = (Ec /(physicoChemical::k.value()*thetaVIntermediate));
            
            label eVibLevel = -1;
            
            if(iMax > SMALL)
            {      
                eVibLevel = cloud_.postCollisionVibrationalEnergyLevel
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

                Ec -= eVibLevel*physicoChemical::k.value()*thetaVIntermediate;
            }                
                            
            // relative translational energy after electronic exchange
            Ec -=  EElistIntermediate[ELevelIntermediate];

            scalar ERot = 0.0;
            
            scalar energyRatio = cloud_.postCollisionRotationalEnergy(rotationalDofIntermediate,ChiBIntermediate);

            ERot = energyRatio*Ec;
        
            Ec -= ERot;
    
            //redistribution finished, test it for dissociation
            
            scalar EcDiss = 0.0;
            label iMax = 0;
            label id = cloud_.constProps(intermediateId_).charDissQuantumLevel()[0];
            
            // calculate if a dissociation is possible
            scalar EVib = eVibLevel*physicoChemical::k.value()*thetaVIntermediate;
            
            EcDiss = Ec + EVib;
            iMax = EcDiss/(physicoChemical::k.value()*thetaVIntermediate);
            
            if((iMax - id) > VSMALL)
            {
                //DISSOCIATION CAN OCCUR
                nReactions_++;
                nReactionsPerTimeStep_++;
                
                if(allowSplitting_)
                {
                    const scalar& heatOfReactionDissociation = heatOfReactionDissociation_*physicoChemical::k.value();
                    
                    translationalEnergy += heatOfReactionDissociation + heatOfReactionRecombinationJoules;
                    
                    translationalEnergy += EEleP + EEleQ + EVibQ + ERotQ;
                    
                    // centre of mass velocity of molecules (pre-split)
                    vector Ucm = (mP*UP + mQ*UQ)/(mP + mQ);
                    
                    mP = cloud_.constProps(productIds_[0]).mass();
                    mQ = cloud_.constProps(productIds_[1]).mass();
                    
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
                    
                    p.typeId() = productIds_[0];
                    p.U() = UP;
                    p.ERot() = 0.0;
                    p.vibLevel().setSize(0,0);
                    p.ELevel() = 0;
                    
                    q.typeId() = productIds_[1];
                    q.U() = UQ;
                    q.ERot() = 0.0;
                    q.vibLevel().setSize(0,0);
                    q.ELevel() = 0;
                }
            }  
        }
    }
}

void  reverseAssociativeIonisation::outputResults(const label& counterIndex)
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
        label nTotReactions = nReactions_;

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
            reduce(nTotReactions, sumOp<label>());
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

        word productMolA = cloud_.typeIdList()[productIds_[0]];
        word productMolB = cloud_.typeIdList()[productIds_[1]];

        if((numberDensities_[0] > 0.0) && (numberDensities_[1] > 0.0))
        {
            scalar reactionRate = 0.0;
            
            reactionRate =
            (
                nTotReactions
                * cloud_.nParticle()
            )/(counterIndex*deltaT*numberDensities_[0]* numberDensities_[1]*volume);
            
            Info<< "Associative ionisation reaction " 
                <<  reactantMolA << " + " << reactantMolB
                <<  " --> " 
                << productMolA << " + " << productMolB
                << ", reaction rate = " << reactionRate
                << endl;
        }
    }
    else
    {
        label nTotReactions = nReactions_;  
        label nReactionsPerTimeStep = nReactionsPerTimeStep_;
        
        if(Pstream::parRun())
        {
            reduce(nTotReactions, sumOp<label>());
            reduce(nReactionsPerTimeStep, sumOp<label>());
        }
        
        if(nTotReactions > VSMALL)
        {
                word reactantMolA = cloud_.typeIdList()[reactantIds_[0]];
                word reactantMolB = cloud_.typeIdList()[reactantIds_[1]];

                word productMolA = cloud_.typeIdList()[productIds_[0]];
                word productMolB = cloud_.typeIdList()[productIds_[1]];
            
                Info<< "Associative ionisation reaction " 
                    <<  reactantMolA << " + " << reactantMolB
                    <<  " --> " 
                    << productMolA << " + " << productMolB
                    << " is active, nReactions this time step = " << nReactionsPerTimeStep << endl;
        } 
    }

    nReactionsPerTimeStep_ = 0.0;
}


const bool& reverseAssociativeIonisation::relax() const
{
    return relax_;
}

}
// End namespace Foam

// ************************************************************************* //
