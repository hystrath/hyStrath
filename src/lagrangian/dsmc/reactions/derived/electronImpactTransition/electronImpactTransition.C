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

#include "electronImpactTransition.H"
#include "addToRunTimeSelectionTable.H"
#include "fvc.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(electronImpactTransition, 0);

addToRunTimeSelectionTable(dsmcReaction, electronImpactTransition, dictionary);



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
electronImpactTransition::electronImpactTransition
(
    Time& t,
    dsmcCloud& cloud,
    const dictionary& dict
)
:
    dsmcReaction(t, cloud, dict),
    propsDict_(dict.subDict(typeName + "Properties")),
    reactantIds_(),
    reactionName_(propsDict_.lookup("reactionName")),   
    relax_(true),
    allowSplitting_(false),
    writeRatesToTerminal_(false),
    volume_(0.0),
    numberDensities_(2, 0.0)
{

}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

electronImpactTransition::~electronImpactTransition()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void electronImpactTransition::initialConfiguration()
{
    setProperties();
}

void electronImpactTransition::setProperties()
{
    // reading in reactants

    const List<word> reactantMolecules (propsDict_.lookup("reactantMolecules"));

    if(reactantMolecules.size() != 2)
    {
        FatalErrorIn("electronImpactTransition::setProperties()")
            << "There should be two or more reactants, instead of " 
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
            FatalErrorIn("electronImpactTransition::setProperties()")
                << "Cannot find type id: " << reactantMolecules[r] << nl 
                << exit(FatalError);
        }
    }    
}

bool electronImpactTransition::tryReactMolecules(const label& typeIdP, const label& typeIdQ)
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

void electronImpactTransition::reaction
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


void electronImpactTransition::reaction
(
    dsmcParcel& p,
    dsmcParcel& q
)
{
    label typeIdP = p.typeId();
    label typeIdQ = q.typeId();

    if(typeIdP == reactantIds_[0] && typeIdQ == reactantIds_[1]) // This produces the correct equilibrium rate A2 + X.
    {

        vector UP = p.U();
        vector UQ = q.U();
//         scalar ERotP = p.ERot();
//         scalar EVibP = p.vibLevel()*cloud_.constProps(typeIdP).thetaV()*physicoChemical::k.value();
        label ELevelP = p.ELevel();//********

        scalar mP = cloud_.constProps(typeIdP).mass();
        scalar mQ = cloud_.constProps(typeIdQ).mass();

        scalar mR = mP*mQ/(mP + mQ);
        scalar cRsqr = magSqr(UP - UQ);
        scalar translationalEnergy = 0.5*mR*cRsqr;

        scalar omegaPQ =
        0.5
        *(
            cloud_.constProps(typeIdP).omega()
            + cloud_.constProps(typeIdQ).omega()
        );

        label jMaxP = cloud_.constProps(typeIdP).numberOfElectronicLevels();//****
        
        List<scalar> EElistP = cloud_.constProps(typeIdP).electronicEnergyList();//***
        
        List<label> gListP = cloud_.constProps(typeIdP).degeneracyList();//***    
        
        label maxLev = 0; // Maximum electronic level energetically possible ****   
        
        label jSelectA = 0; // selected intermediate integer electronic level (0 to jSelect).
        label jSelectB = 0; // selected intermediate integer electronic level (0 to jSelect).
        label jSelect = 0; // Minimum of jSelectA and jSelectB. 
        
        scalar g = 0.0; // Distribution function of Eq. 3.1.6 of Liechty thesis.
        scalar gMax = 0.0; // Maximum value of Eq. 3.1.6     
        
        label jDash = 0; // random integer between 0 and jSelect.
        scalar EJ = 0.0; // maximum possible electronic energy level within list based on EcP.
        label gJ = 0; // maximum possible degeneracy level within list.
        scalar denomMax = 0.0; // maximum denominator value in Liechty pdf (see below).
        scalar func = 0.0; // distribution function Eq. 3.1.8 of Liechty thesis.

        scalar EcPP = 0.0;
        
        // calculate if electronImpactTransition of species P is possible
        EcPP = translationalEnergy + EElistP[ELevelP]; //*****

        for (label ii = 0; ii < jMaxP; ++ii)
        {  
            if ((EElistP[ii]) > EcPP) break;
            maxLev = ii;
            jSelectA = ii;
            
            //Eq. 3.1.6 of Liechty thesis.	  
            g = gListP[ii]*pow((EcPP - EElistP[ii]),(1.5 - omegaPQ));
            
            if ( ii == 0 || gMax < g )
            {
                gMax = g;
                jSelectB = ii;
            }
        }
        
        jSelect = jSelectA;
        if (jSelectB < jSelectA) jSelect = jSelectB; 

        EJ = EElistP[jSelect]; //Max. poss energy in list : list goes from [0] to [jMaxP-1]
        gJ = gListP[jSelect]; //Max. poss degeneracy in list : list goes from [0] to [jMaxP-1]
        denomMax = gJ*pow((EcPP - EJ), (1.5 - omegaPQ)); // Max. in denominator of Liechty pdf for post-collision pdf.
    
        do // acceptance - rejection based on Eq. 3.1.8 of Liechty thesis.	
        {
            jDash = cloud_.rndGen().integer(0,maxLev);
            func = gListP[jDash]*pow((EcPP - EElistP[jDash]), (1.5 - omegaPQ))/denomMax;

        } while( !(func > cloud_.rndGen().sample01<scalar>()));

        if (ELevelP == 1 && jDash == 2)
        {
            nTotReactions_++;
        }
    }
  
    if(typeIdP == reactantIds_[1] && typeIdQ == reactantIds_[0]) // This produces the correct equilibrium rate X + A2.
    {

        vector UP = p.U();
        vector UQ = q.U();
//         scalar ERotQ = q.ERot();
//         scalar EVibQ = q.vibLevel()*cloud_.constProps(typeIdQ).thetaV()*physicoChemical::k.value();
        label ELevelQ = q.ELevel();//********

        scalar mP = cloud_.constProps(typeIdP).mass();
        scalar mQ = cloud_.constProps(typeIdQ).mass();

        scalar mR = mP*mQ/(mP + mQ);
        scalar cRsqr = magSqr(UP - UQ);
        scalar translationalEnergy = 0.5*mR*cRsqr;

        scalar omegaPQ =
        0.5
        *(
            cloud_.constProps(typeIdP).omega()
            + cloud_.constProps(typeIdQ).omega()
        );

        label jMaxQ = cloud_.constProps(typeIdQ).numberOfElectronicLevels();//****
        
        List<scalar> EElistQ = cloud_.constProps(typeIdQ).electronicEnergyList();//***
        
        List<label> gListQ = cloud_.constProps(typeIdQ).degeneracyList();//***    
        
        label maxLev = 0; // Maximum electronic level energetically possible ****   
        
        label jSelectA = 0; // selected intermediate integer electronic level (0 to jSelect).
        label jSelectB = 0; // selected intermediate integer electronic level (0 to jSelect).
        label jSelect = 0; // Minimum of jSelectA and jSelectB. 
        
        scalar g = 0.0; // Distribution function of Eq. 3.1.6 of Liechty thesis.
        scalar gMax = 0.0; // Maximum value of Eq. 3.1.6     
        
        label jDash = 0; // random integer between 0 and jSelect.
        scalar EJ = 0.0; // maximum possible electronic energy level within list based on EcP.
        label gJ = 0; // maximum possible degeneracy level within list.
        scalar denomMax = 0.0; // maximum denominator value in Liechty pdf (see below).
        scalar func = 0.0; // distribution function Eq. 3.1.8 of Liechty thesis.      
        scalar EcPQ = 0.0;
            
        // calculate if electronImpactTransition of species P is possible
        EcPQ = translationalEnergy + EElistQ[ELevelQ]; //*****

        for (label ii = 0; ii < jMaxQ; ++ii)
        {
            
            if ((EElistQ[ii]) > EcPQ) break;
            
            maxLev = ii;
            jSelectA = ii;
            
            //Eq. 3.1.6 of Liechty thesis.	  
            g = gListQ[ii]*pow((EcPQ - EElistQ[ii]),(1.5 - omegaPQ));

            if ( ii == 0 || gMax < g )
            {
                gMax = g;
                jSelectB = ii;
            }
        }
        
        jSelect = jSelectA;
        if (jSelectB < jSelectA) jSelect = jSelectB;

        EJ = EElistQ[jSelect]; //Max. poss energy in list : list goes from [0] to [jMaxP-1]
        gJ = gListQ[jSelect]; //Max. poss degeneracy in list : list goes from [0] to [jMaxP-1]
        denomMax = gJ*pow((EcPQ - EJ), (1.5 - omegaPQ)); // Max. in denominator of Liechty pdf for post-collision pdf.
     
        do // acceptance - rejection based on Eq. 3.1.8 of Liechty thesis.	
        {
            jDash = cloud_.rndGen().integer(0,maxLev);
            func = gListQ[jDash]*pow((EcPQ - EElistQ[jDash]), (1.5 - omegaPQ))/denomMax;

        } while( !(func > cloud_.rndGen().sample01<scalar>()));

        if (ELevelQ == 1 && jDash == 2)
        {
            nTotReactions_++;
        }
     }
}

void  electronImpactTransition::outputResults(const label& counterIndex)
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
                    if(id == 0)
                    {
                        if(p->ELevel() == 1)
                        {
                            mols[id]++;
                        }
                    }
                    else
                    {
                        mols[id]++;
                    }
                }
            }

            volume_ += mesh_.cellVolumes()[c];
        }
        
        scalar volume = volume_;
        label nTotReactions = nTotReactions_;

        //- Parallel communication
        if(Pstream::parRun())
        {
            reduce(volume, sumOp<scalar>());
            reduce(mols[0], sumOp<label>());
            reduce(mols[1], sumOp<label>());
            reduce(nTotReactions, sumOp<label>());
        }

        numberDensities_[0] = (mols[0]*cloud().nParticle())/volume;
        numberDensities_[1] = (mols[1]*cloud().nParticle())/volume;

        const scalar& deltaT = mesh_.time().deltaT().value();

        word reactantMolA = cloud_.typeIdList()[reactantIds_[0]];
        word reactantMolB = cloud_.typeIdList()[reactantIds_[1]];

        if((numberDensities_[0] > 0.0) && (numberDensities_[1] > 0.0))
        { 

            scalar reactionRate = 0.0;

            reactionRate =
            (
            nTotReactions
            * cloud_.nParticle()
            )/(counterIndex*deltaT*numberDensities_[0]* numberDensities_[1]*volume);

/*            Info<< "cloud_.nParticle() = "<<  cloud_.nParticle() << endl;
        Info<< "counterIndex = "<<  counterIndex << endl;
        Info<< "deltaT = "<<  deltaT << endl;
        Info<< " numberDensities_[0] = "<< numberDensities_[0]  << endl;
        Info<< " numberDensities_[1] = "<< numberDensities_[1]  << endl;
        Info<< " volume = "<< volume  << endl;*/	    
        
            Info<< "electronImpactTransition "
                <<  reactantMolA << " + " << reactantMolB 
                <<  " --> "
                << reactantMolA << " + " << reactantMolB 
                << ", reaction rate = " << reactionRate
                << endl;
        }
    }
    else
    {
        label nTotReactions = nTotReactions_;   
        
        if(Pstream::parRun())
        {
            reduce(nTotReactions, sumOp<label>());
        }
    
        if(nTotReactions > VSMALL)
        {
                word reactantMolA = cloud_.typeIdList()[reactantIds_[0]];
                word reactantMolB = cloud_.typeIdList()[reactantIds_[1]];
            
                Info<< "electronImpactTransition reaction "
                <<  reactantMolA << " + " << reactantMolB 
                <<  " --> "
                << reactantMolA << " + " << reactantMolB  
                    << " is active." << endl;
        }       
    }

    nReactionsPerTimeStep_ = 0.0;

}


const bool& electronImpactTransition::relax() const
{
    return relax_;
}

}
// End namespace Foam

// ************************************************************************* //
