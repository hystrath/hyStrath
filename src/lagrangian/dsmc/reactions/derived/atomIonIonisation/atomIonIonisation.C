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

#include "atomIonIonisation.H"
#include "addToRunTimeSelectionTable.H"
#include "fvc.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(atomIonIonisation, 0);

addToRunTimeSelectionTable(dsmcReaction, atomIonIonisation, dictionary);



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
atomIonIonisation::atomIonIonisation
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
    heatOfReactionIon_(readScalar(propsDict_.lookup("heatOfReactionIonisation"))),
    nTotIonisationReactions_(0),
    nIonisationReactionsPerTimeStep_(0),
    relax_(true),
    allowSplitting_(false),
    writeRatesToTerminal_(false),
    volume_(0.0),
    numberDensities_(2, 0.0)
{

}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

atomIonIonisation::~atomIonIonisation()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void atomIonIonisation::initialConfiguration()
{
    setProperties();
}

void atomIonIonisation::setProperties()
{
    // reading in reactants

    const List<word> reactantMolecules (propsDict_.lookup("reactants"));

    if(reactantMolecules.size() != 2)
    {
        FatalErrorIn("atomIonIonisation::setProperties()")
            << "There should be two or more reactants, instead of " 
            << reactantMolecules.size() << nl 
            << exit(FatalError);
    }
    
    if(reactantMolecules[0] == reactantMolecules[1])
    {
        FatalErrorIn("atomIonIonisation::setProperties()")
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
            FatalErrorIn("atomIonIonisation::setProperties()")
                << "Cannot find type id: " << reactantMolecules[r] << nl 
                << exit(FatalError);
        }
    }
    
    // check that reactant one is an 'ATOM' 

    const label& rDof1 = cloud_.constProps(reactantIds_[0]).rotationalDegreesOfFreedom();

    if(rDof1 > VSMALL)
    {
        FatalErrorIn("atomIonIonisation::setProperties()")
            << "First reactant must be an atom (not a molecule or an electron): " << reactantMolecules[0] 
            << nl 
            << exit(FatalError);
    }
    
    // check that reactant two is a charged 'MOLECULE'

    const label& rDof2 = cloud_.constProps(reactantIds_[1]).rotationalDegreesOfFreedom();

    if(rDof2 < VSMALL)
    {
        FatalErrorIn("atomIonIonisation::setProperties()")
            << "Second reactant must be a charged molecule: " << reactantMolecules[1] 
            << nl 
            << exit(FatalError);
    }
    
    const label& charge2 = cloud_.constProps(reactantIds_[1]).charge();

    if(charge2 != 1)
    {
        FatalErrorIn("atomIonIonisation::setProperties()")
            << "Second reactant must be a charged molecule: " << reactantMolecules[1] 
            << nl 
            << exit(FatalError);
    }
    
    // check that reactant two only has a single vibrational degree of freedom 

    const label& vDof = cloud_.constProps(reactantIds_[1]).vibrationalDegreesOfFreedom();

    if(vDof > 1)
    {
        FatalErrorIn("atomIonIonisation::setProperties()")
            << "Reactions are currently only implemented for monatomic and diatomic species"
            << " This is a polyatomic:" << reactantMolecules[1] 
            << nl 
            << exit(FatalError);
    }
    
//     if(chargedAtom_)
//     {
        // reading in ionisation products

        List<word> productMoleculesIonisation (propsDict_.lookup("productsOfIonisedAtom"));

        if(productMoleculesIonisation.size() != 2)
        {
            FatalErrorIn("atomIonIonisation::setProperties()")
                << "Number of ionisation products is " << productMoleculesIonisation.size() <<
                ", should be two."
                << exit(FatalError);
        }
        

        productIdsIon_.setSize(productMoleculesIonisation.size());

        forAll(productMoleculesIonisation, r)
        {
            if(productIdsIon_.size() != 2)
            {
                FatalErrorIn("atomIonIonisation::setProperties()")
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
                    FatalErrorIn("atomIonIonisation::setProperties()")
                        << "Cannot find type id: " << productMoleculesIonisation[r] << nl 
                        << exit(FatalError);
                }
            }
            
            // check that product one is a 'ATOM', not an 'MOLECULE'

            const scalar& rDof3 = cloud_.constProps(productIdsIon_[0]).rotationalDegreesOfFreedom();

            if(rDof3 > 1)
            {
                FatalErrorIn("atomIonIonisation::setProperties()")
                    << "First product must be an atom (not an atom): " << productMoleculesIonisation[0] 
                    << nl 
                    << exit(FatalError);
            }

            // check that product two is an 'ELECTRON'

            const label& charge = cloud_.constProps(productIdsIon_[1]).charge();

            if(charge != -1)
            {
                FatalErrorIn("atomIonIonisation::setProperties()")
                    << "Second product must be an electron: " << productMoleculesIonisation[1] 
                    << nl 
                    << exit(FatalError);
            }
        }
//     }
}

bool atomIonIonisation::tryReactMolecules(const label& typeIdP, const label& typeIdQ)
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

void atomIonIonisation::reaction
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


void atomIonIonisation::reaction
(
    dsmcParcel& p,
    dsmcParcel& q
)
{
    label typeIdP = p.typeId();
    label typeIdQ = q.typeId();

    if( typeIdP == reactantIds_[0] && typeIdQ == reactantIds_[1])
    {
        relax_ = true;
        
        vector UP = p.U();
        vector UQ = q.U();
        scalar ERotP = p.ERot();
        scalar ERotQ = q.ERot();
        scalar EVibQ = 0;
        if(cloud_.constProps(typeIdQ).vibrationalDegreesOfFreedom() > 0)
        {        
            EVibQ = q.vibLevel()[0]
                        *cloud_.constProps(typeIdQ).thetaV()[0]
                        *physicoChemical::k.value();
        }
        scalar EEleP = cloud_.constProps(typeIdP).electronicEnergyList()[p.ELevel()];
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

        scalar heatOfReactionJoulesIon = heatOfReactionIon_*physicoChemical::k.value();

        
        scalar omegaPQ =
            0.5
            *(
                    cloud_.constProps(typeIdP).omega()
                + cloud_.constProps(typeIdQ).omega()
            );
            
        scalar ChiB = 2.5 - omegaPQ;
        
        scalar EcPPIon= 0.0;
        scalar ionisationEnergy = cloud_.constProps(typeIdP).ionisationTemperature()*physicoChemical::k.value();
        
        // calculate if an ionisation of species P is possible
        EcPPIon = translationalEnergy + EEleP;
        
        if((EcPPIon - ionisationEnergy) > VSMALL)
        {
            nTotIonisationReactions_++;
            nIonisationReactionsPerTimeStep_++;
            
            if(allowSplitting_)
            {
                relax_ = false;
                
                scalar thetaVQ = 0;
                scalar thetaDQ = 0;
                scalar ZrefQ = 0;
                scalar refTempZvQ = 0;
                
                if(cloud_.constProps(typeIdQ).vibrationalDegreesOfFreedom() > 0)
                {        
                    thetaVQ = cloud_.constProps(typeIdQ).thetaV()[0];
                    thetaDQ = cloud_.constProps(typeIdQ).thetaD()[0];
                    ZrefQ = cloud_.constProps(typeIdQ).Zref()[0];
                    refTempZvQ = cloud_.constProps(typeIdQ).TrefZv()[0];
                }
                scalar jMaxQ = cloud_.constProps(typeIdQ).numberOfElectronicLevels()-1;
                scalar rotationalDofQ = cloud_.constProps(typeIdQ).rotationalDegreesOfFreedom();
                
                                
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
                
                label vibLevelQ = 0;
                
                if(rotationalDofQ > 0)
                {
                    translationalEnergy += EVibQ;
                    
                    label iMax = (translationalEnergy / (physicoChemical::k.value()*thetaVQ));
                    
                    vibLevelQ = cloud_.postCollisionVibrationalEnergyLevel
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
                }
                
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
                
                scalar translationalEnergy2 = ERotP;

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
                q.ERot() = ERotQ;
                q.vibLevel()[0] = vibLevelQ;

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
                    -1,
                    classificationP,
                    vibLevel
                );
            }
        }
    }
    
    if(typeIdP == reactantIds_[1] && typeIdQ == reactantIds_[0]) // This produces the correct equilibrium rate A2 + X.
    {   
        relax_ = true;
        
        vector UP = p.U();
        vector UQ = q.U();
        scalar ERotP = p.ERot();
        scalar ERotQ = q.ERot();
        scalar EVibP = 0.0;
        if(cloud_.constProps(typeIdP).vibrationalDegreesOfFreedom() > 0)
        {        
            EVibP = p.vibLevel()[0]
                        *cloud_.constProps(typeIdP).thetaV()[0]
                        *physicoChemical::k.value();
        }
        scalar EEleP = cloud_.constProps(typeIdP).electronicEnergyList()[p.ELevel()];
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

        scalar heatOfReactionJoulesIon = heatOfReactionIon_*physicoChemical::k.value();
       
        
        scalar omegaPQ =
            0.5
            *(
                    cloud_.constProps(typeIdP).omega()
                + cloud_.constProps(typeIdQ).omega()
            );
            
        scalar ChiB = 2.5 - omegaPQ;
        
        scalar EcPPIon= 0.0;
        
        scalar ionisationEnergy = cloud_.constProps(typeIdQ).ionisationTemperature()*physicoChemical::k.value();
        
        // calculate if an ionisation of species Q is possible
        EcPPIon = translationalEnergy + EEleQ;
        
        if((EcPPIon - ionisationEnergy) > VSMALL)
        {
            nTotIonisationReactions_++;
            nIonisationReactionsPerTimeStep_++;
            
            if(allowSplitting_)
            {
                relax_ = false;
                
                scalar thetaVP = 0;
                scalar thetaDP = 0;
                scalar ZrefP = 0;
                scalar refTempZvP = 0;
                
                if(cloud_.constProps(typeIdP).vibrationalDegreesOfFreedom() > 0)
                {        
                    thetaVP = cloud_.constProps(typeIdP).thetaV()[0];
                    thetaDP = cloud_.constProps(typeIdP).thetaD()[0];
                    ZrefP = cloud_.constProps(typeIdP).Zref()[0];
                    refTempZvP = cloud_.constProps(typeIdP).TrefZv()[0];
                }
                scalar jMaxP = cloud_.constProps(typeIdP).numberOfElectronicLevels()-1;
                scalar rotationalDofP = cloud_.constProps(typeIdP).rotationalDegreesOfFreedom();
               
                                
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
                
                label vibLevelP = 0;
                
                if(rotationalDofP > 0)
                {
                    translationalEnergy += EVibP;
                    
                    label iMax = (translationalEnergy / (physicoChemical::k.value()*thetaVP));
                    
                    vibLevelP = cloud_.postCollisionVibrationalEnergyLevel
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
                }
                
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
                
                scalar translationalEnergy2 = ERotQ;

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
                p.ERot() = ERotP;
                p.vibLevel()[0] = vibLevelP;

                // Atom Q will ionise
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
    }
}

void  atomIonIonisation::outputResults(const label& counterIndex)
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
        label nTotReactionsIonisation = nTotIonisationReactions_;

        //- Parallel communication
        if(Pstream::parRun())
        {
            reduce(volume, sumOp<scalar>());
            reduce(mols[0], sumOp<label>());
            reduce(mols[1], sumOp<label>());
            reduce(nTotReactionsIonisation, sumOp<label>());
        }

        numberDensities_[0] = (mols[0]*cloud().nParticle())/volume;
        numberDensities_[1] = (mols[1]*cloud().nParticle())/volume;

        const scalar& deltaT = mesh_.time().deltaT().value();

        word reactantMolA = cloud_.typeIdList()[reactantIds_[0]];
        word reactantMolB = cloud_.typeIdList()[reactantIds_[1]];
        
        word productMolC = cloud_.typeIdList()[productIdsIon_[0]];
        word productMolD = cloud_.typeIdList()[productIdsIon_[1]];

        if((numberDensities_[0] > 0.0) && (numberDensities_[1] > 0.0))
        {   
            scalar reactionRateIonisation = 0.0;

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
        }
    }
    else
    {
        label nTotReactionsIonisation = nTotIonisationReactions_;   
        label nIonisationReactionsPerTimeStep = nIonisationReactionsPerTimeStep_;
        
        if(Pstream::parRun())
        {
            reduce(nTotReactionsIonisation, sumOp<label>());
            reduce(nIonisationReactionsPerTimeStep, sumOp<label>());
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
    }

    nIonisationReactionsPerTimeStep_ = 0.0;
}


const bool& atomIonIonisation::relax() const
{
    return relax_;
}

}
// End namespace Foam

// ************************************************************************* //
