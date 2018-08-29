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

#include "dissociationTypeIDissimilarSpecies.H"
#include "addToRunTimeSelectionTable.H"
#include "fvc.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(dissociationTypeIDissimilarSpecies, 0);

addToRunTimeSelectionTable(dsmcReaction, dissociationTypeIDissimilarSpecies, dictionary);



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
dissociationTypeIDissimilarSpecies::dissociationTypeIDissimilarSpecies
(
    Time& t,
    dsmcCloud& cloud,
    const dictionary& dict
)
:
    dsmcReaction(t, cloud, dict),
    propsDict_(dict.subDict(typeName + "Properties")),
    reactantIds_(),
    productsToDiss_(),
    reactionName_(propsDict_.lookup("reactionName")),
    heatOfReactionAB_(readScalar(propsDict_.lookup("heatOfReactionAB"))),
    heatOfReactionCD_(readScalar(propsDict_.lookup("heatOfReactionCD"))),
    nABReactions_(0),
    nCDReactions_(0),
    relax_(true),
    allowSplitting_(false),
    writeRatesToTerminal_(false),
    volume_(0.0),
    numberDensities_(2, 0.0)
{

}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

dissociationTypeIDissimilarSpecies::~dissociationTypeIDissimilarSpecies()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void dissociationTypeIDissimilarSpecies::initialConfiguration()
{
    setProperties();
}

void dissociationTypeIDissimilarSpecies::setProperties()
{
    // reading in reactants

    const List<word> reactantMolecules (propsDict_.lookup("reactantMoleculesToDissociate"));

    if(reactantMolecules.size() != 2)
    {
        FatalErrorIn("dissociationTypeIDissimilarSpecies::setProperties()")
            << "There should be two or more reactants, instead of " 
            << reactantMolecules.size() << nl 
            << exit(FatalError);
    }
    
    if(reactantMolecules[0] == reactantMolecules[1])
    {
        FatalErrorIn("dissociationTypeIDissimilarSpecies::setProperties()")
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
            FatalErrorIn("dissociationTypeIDissimilarSpecies::setProperties()")
                << "Cannot find type id: " << reactantMolecules[r] << nl 
                << exit(FatalError);
        }

        // check that reactants are 'MOLECULES' (not 'ATOMS') 

        const scalar& rDof = cloud_.constProps(reactantIds_[r]).rotationalDegreesOfFreedom();
    
        if(rDof < 1)
        {
            FatalErrorIn("dissociationTypeIDissimilarSpecies::setProperties()")
                << "Reactant must be a molecule (not an atom): " << reactantMolecules[r] 
                << nl 
                << exit(FatalError);
        }
    }

    // reading in products

    List< List<word> > productMolecules (propsDict_.lookup("productsOfDissociatedMolecule"));

    if(productMolecules.size() != reactantIds_.size())
    {
        FatalErrorIn("dissociationTypeIDissimilarSpecies::setProperties()")
            << "number of reactant molecules to be dissociated = " << reactantIds_.size()
            << " is not the same as the number of products = " << productMolecules.size()
            << exit(FatalError);
    }
    

    productsToDiss_.setSize(productMolecules.size());

    forAll(productMolecules, r)
    {
        const List<word>& productsForDiss = productMolecules[r];

        if(productsForDiss.size() != 2)
        {
            FatalErrorIn("dissociationTypeIDissimilarSpecies::setProperties()")
                << "There should be two products (for the dissociating molecule "
                << reactantMolecules[r] << "), instead of " 
                << productsForDiss.size() << nl 
                << exit(FatalError);
        }
    
        productsToDiss_[r].setSize(productsForDiss.size(), -1);
    
        forAll(productsToDiss_[r], p)
        {
            productsToDiss_[r][p] = findIndex(cloud_.typeIdList(), productsForDiss[p]);
        
            if(productsToDiss_[r][p] == -1)
            {
                FatalErrorIn("dissociationTypeIDissimilarSpecies::setProperties()")
                    << "Cannot find type id: " << productsForDiss[p] << nl 
                    << exit(FatalError);
            }
        }
    }
}



bool dissociationTypeIDissimilarSpecies::tryReactMolecules(const label& typeIdP, const label& typeIdQ)
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

void dissociationTypeIDissimilarSpecies::reaction
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
//£££££££££££££££££££££££££££££££££££££££££££££££££££££££££££££££££££££££££££££££££££££££££££
void dissociationTypeIDissimilarSpecies::reaction
(
    dsmcParcel& p,
    dsmcParcel& q
)

{
    label typeIdP = p.typeId();
    label typeIdQ = q.typeId();

    if(typeIdP == reactantIds_[0] && typeIdQ == reactantIds_[1])
    {
        vector UP = p.U();
        vector UQ = q.U();
        scalar ERotP = p.ERot();
        scalar ERotQ = q.ERot();
        scalar EVibP = p.vibLevel()*cloud_.constProps(typeIdP).thetaV()*physicoChemical::k.value();
        scalar EVibQ = q.vibLevel()*cloud_.constProps(typeIdQ).thetaV()*physicoChemical::k.value();

        scalar mP = cloud_.constProps(typeIdP).mass();
        scalar mQ = cloud_.constProps(typeIdQ).mass();

        scalar mR = mP*mQ/(mP + mQ);
        scalar cRsqr = magSqr(UP - UQ);
        scalar translationalEnergy = 0.5*mR*cRsqr;
        
        scalar thetaVP = cloud_.constProps(typeIdP).thetaV();
        scalar thetaVQ = cloud_.constProps(typeIdQ).thetaV();
        
        scalar thetaDP = cloud_.constProps(typeIdP).thetaD();
        scalar thetaDQ = cloud_.constProps(typeIdQ).thetaD();
        
        scalar ZrefP = cloud_.constProps(typeIdP).Zref();
        scalar ZrefQ = cloud_.constProps(typeIdQ).Zref();
        
        scalar refTempZvP = cloud_.constProps(typeIdP).TrefZv();
        scalar refTempZvQ = cloud_.constProps(typeIdQ).TrefZv();
        
        scalar rotationalDofP = cloud_.constProps(typeIdP).rotationalDegreesOfFreedom();
        scalar rotationalDofQ = cloud_.constProps(typeIdQ).rotationalDegreesOfFreedom();

        scalar omegaPQ =
        0.5
        *(
            cloud_.constProps(typeIdP).omega()
            + cloud_.constProps(typeIdQ).omega()
        );
        
        scalar ChiB = 2.5 - omegaPQ;

        bool dissoP = false;
        bool dissoQ = false;

        scalar EcPQ = 0.0;
        label idP = cloud_.constProps(typeIdP).charDissQuantumLevel();
        label imaxP = 0;

        // calculate if a dissociation of species P is possible
        EcPQ = translationalEnergy + EVibP;

        imaxP = EcPQ/(physicoChemical::k.value()*cloud_.constProps(typeIdP).thetaV());
        
        if(imaxP-idP > 0)
        {
            dissoP = true; 
        }

        scalar EcQP = 0.0;
        label idQ = cloud_.constProps(typeIdQ).charDissQuantumLevel();
        label imaxQ = 0;

        // calculate if a dissociation of species Q is possible
        EcQP = translationalEnergy + EVibQ;
        
        if(dissoP == true)
        {
            EcQP += (heatOfReactionAB_*physicoChemical::k.value());
        }

        imaxQ = EcQP/(physicoChemical::k.value()*cloud_.constProps(typeIdQ).thetaV());

        
        if(imaxQ-idQ > 0) 
        {
            dissoQ = true;
        }
        
        relax_ = true;

        if(dissoP == true && dissoQ == false)
        {
            nABReactions_++;

            relax_ = false; 
    
            if(allowSplitting_)
            {
                scalar heatOfReactionJoules = heatOfReactionAB_*physicoChemical::k.value();
                
                translationalEnergy = translationalEnergy + heatOfReactionJoules + EVibP;
                
//                 if (0.2 > cloud_.rndGen().sample01<scalar>())
//                 {
//                     scalar EcP = translationalEnergy + ERotP;
//                     
//                     scalar energyRatio = 0.0;
//                     
//                     if(rotationalDofP == 2.0)
//                     {
//                         energyRatio = 1.0 - pow(cloud_.rndGen().sample01<scalar>(),(1.0/ChiB));
//                     }
//                     else
//                     {
//                         scalar ChiA = 0.5*rotationalDofP;
//                         
//                         energyRatio = cloud_.energyRatio(ChiA, ChiB);
//                     }
// 
//                     ERotP = energyRatio*EcP;
//                 
//                     translationalEnergy = EcP - ERotP;
//                 }
                
                scalar EcQ = translationalEnergy + EVibQ;
                label iMaxQ = (EcQ /(physicoChemical::k.value()*thetaVQ));

                if(iMaxQ > SMALL)
                {
                    // - Bird equations 5.42 and 11.34 gave this denominator
                    scalar TCollQ = (iMaxQ*thetaVQ) / (3.5 - omegaPQ); 
                    
                    scalar pow1 = pow((thetaDQ/TCollQ),0.333) - 1.0;

                    scalar pow2 = pow ((thetaDQ/refTempZvQ),0.333) -1.0;
                    
                    // - vibrational collision number (equation 2, Bird 2010)
                    scalar ZvQ1 = pow((thetaDQ/TCollQ),omegaPQ); 
                    
                    scalar ZvQ2 = pow(ZrefQ*(pow((thetaDQ/refTempZvQ),(-1.0*omegaPQ))),(pow1/pow2));
                    
                    scalar ZvQ = ZvQ1*ZvQ2;
//                     scalar ZvQ = 50.0;
                        
                    scalar inverseVibrationalCollisionNumberQ = 1.0/(5.0*ZvQ);
                
                    if(inverseVibrationalCollisionNumberQ > cloud_.rndGen().sample01<scalar>())
                    {
                        label iDashQ = 0; // post-collision quantum number
                        scalar func = 0.0;

                        do // acceptance - rejection 
                        {
                            iDashQ = cloud_.rndGen().integer(0,iMaxQ);
                            EVibQ = iDashQ*physicoChemical::k.value()*thetaVQ;
                            func = pow((1.0 - (EVibQ / EcQ)),(1.5 - omegaPQ));
                    
                        } while( !(func > cloud_.rndGen().sample01<scalar>()) );
                
                        translationalEnergy = EcQ - EVibQ;
                    }
                }
                
                if (0.2 > cloud_.rndGen().sample01<scalar>())
                {
                    scalar EcQ = translationalEnergy + ERotQ;
                    
                    scalar energyRatio = 0.0;
                    
                    if(rotationalDofQ == 2.0)
                    {
                        energyRatio = 1.0 - pow(cloud_.rndGen().sample01<scalar>(),(1.0/ChiB));
                    }
                    else
                    {
                        scalar ChiA = 0.5*rotationalDofQ;
                        
                        energyRatio = cloud_.energyRatio(ChiA, ChiB);
                    }

                    ERotQ = energyRatio*EcQ;
                
                    translationalEnergy = EcQ - ERotQ;
                }
                
                scalar relVelNonDissoMol = sqrt((2.0*translationalEnergy)/mR);

                //center of mass velocity of all reactant particles
                vector Ucm = (mP*UP + mQ*UQ)/(mP + mQ);

                // Variable Hard Sphere collision part for collision of molecules
    
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
    
                UP = Ucm + (postCollisionRelU*mQ/(mP + mQ)); // UP is used as Ucm for atomic split.

                UQ = Ucm - (postCollisionRelU*mP/(mP + mQ)); // Q is the NON-DISSOCIATING molecule.

                // Q remains NON-DISSOCIATED (but internal and translational energy modified).

                label p1 = 0;
                label p2 = 1;

                const label& typeId1 = productsToDiss_[0][p1]; //note [1] corresponds to the second set of products in chemReactDict
                const label& typeId2 = productsToDiss_[0][p2];

                //Mass of Product one and two
                scalar mP1 = cloud_.constProps(typeId1).mass();
                scalar mP2 = cloud_.constProps(typeId2).mass();

                // Q remains NON-DISSOCIATED (but internal and translational energy modified).
                q.vibLevel() = EVibQ/(cloud_.constProps(typeIdQ).thetaV()*physicoChemical::k.value());
                q.ERot() = ERotQ;
                q.U() = UQ;

                scalar mRatoms = mP1*mP2/(mP1 + mP2);

                vector UcmAtomsP = UP;

//                 scalar cRatoms = sqrt(2.0*Ecplx/mRatoms);
                
                scalar cRatoms = sqrt(2.0*ERotP/mRatoms);

                // Variable Hard Sphere collision part for atoms
            
                cosTheta = 2.0*cloud_.rndGen().sample01<scalar>() - 1.0;
            
                sinTheta = sqrt(1.0 - cosTheta*cosTheta);
            
                phi = twoPi*cloud_.rndGen().sample01<scalar>();
            
                vector postCollisionRelU2 = cRatoms
                *vector
                    (
                        cosTheta,
                        sinTheta*cos(phi),
                        sinTheta*sin(phi)
                    );

                vector uP1 = UcmAtomsP + (postCollisionRelU2*mP2/(mP1 + mP2));
                vector uP2 = UcmAtomsP - (postCollisionRelU2*mP1/(mP1 + mP2));

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
                p.vibLevel() = 0;
                p.ERot() = 0.0;
                
                label classificationP = p.classification();
                scalar RWF = p.RWF();
                
                // insert new product 2
                cloud_.addNewParcel
                (
                    position,
                    uP2,
                    RWF,
                    0.0,
                    0,
                    0,
                    cell,
                    tetFace,
                    tetPt,
                    typeId2,
                    -1,
                    classificationP
                );
            }
        }

        if(dissoQ == true && dissoP == false)
        {
            nCDReactions_++;

            relax_ = false;
            
            if(allowSplitting_)
            {
                scalar heatOfReactionJoules = heatOfReactionCD_*physicoChemical::k.value();
                
                translationalEnergy = translationalEnergy + heatOfReactionJoules + EVibQ;   
                
                scalar EcP = translationalEnergy + EVibP;
                label iMaxP = (EcP /(physicoChemical::k.value()*thetaVP));

                if(iMaxP > SMALL)
                {
                    // - Bird equations 5.42 and 11.34 gave this denominator
                    scalar TCollP = (iMaxP*thetaVP) / (3.5 - omegaPQ); 
                    
                    scalar pow1 = pow((thetaDP/TCollP),0.333) - 1.0;

                    scalar pow2 = pow ((thetaDP/refTempZvP),0.333) -1.0;
                    
                    // - vibrational collision number (equation 2, Bird 2010)
                    scalar ZvP1 = pow((thetaDP/TCollP),omegaPQ); 
                    
                    scalar ZvP2 = pow(ZrefP*(pow((thetaDP/refTempZvP),(-1.0*omegaPQ))),(pow1/pow2));
                    
                    scalar ZvP = ZvP1*ZvP2;
//                     scalar ZvP = 50.0;
                        
                    scalar inverseVibrationalCollisionNumberQ = 1.0/(5.0*ZvP);
                
                    if(inverseVibrationalCollisionNumberQ > cloud_.rndGen().sample01<scalar>())
                    {
                        label iDashP = 0; // post-collision quantum number
                        scalar func = 0.0;

                        do // acceptance - rejection 
                        {
                            iDashP = cloud_.rndGen().integer(0,iMaxP);
                            EVibP = iDashP*physicoChemical::k.value()*thetaVP;
                            func = pow((1.0 - (EVibP / EcP)),(1.5 - omegaPQ));
                    
                        } while( !(func > cloud_.rndGen().sample01<scalar>()) );
                
                        translationalEnergy = EcP - EVibP;
                    }
                }
                
                if (0.2 > cloud_.rndGen().sample01<scalar>())
                {
                    scalar EcP = translationalEnergy + ERotP;
                    
                    scalar energyRatio = 0.0;
                    
                    if(rotationalDofP == 2.0)
                    {
                        energyRatio = 1.0 - pow(cloud_.rndGen().sample01<scalar>(),(1.0/ChiB));
                    }
                    else
                    {
                        scalar ChiA = 0.5*rotationalDofP;
                        
                        energyRatio = cloud_.energyRatio(ChiA, ChiB);
                    }

                    ERotP = energyRatio*EcP;
                
                    translationalEnergy = EcP - ERotP;
                }
                
                //Molecule Q is dissociating, so no vibrational relaxation to it
                
//                 if (0.2 > cloud_.rndGen().sample01<scalar>())
//                 {
//                     scalar EcQ = translationalEnergy + ERotQ;
//                     
//                     scalar energyRatio = 0.0;
//                     
//                     if(rotationalDofQ == 2.0)
//                     {
//                         energyRatio = 1.0 - pow(cloud_.rndGen().sample01<scalar>(),(1.0/ChiB));
//                     }
//                     else
//                     {
//                         scalar ChiA = 0.5*rotationalDofQ;
//                         
//                         energyRatio = cloud_.energyRatio(ChiA, ChiB);
//                     }
// 
//                     ERotQ = energyRatio*EcQ;
//                 
//                     translationalEnergy = EcQ - ERotQ;
//                 }
                
                scalar relVelNonDissoMol = sqrt((2.0*translationalEnergy)/mR);

                //center of mass velocity of all reactant particles
                vector Ucm = (mP*UP + mQ*UQ)/(mP + mQ);

                // Variable Hard Sphere collision part for collision of molecules
        
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
        
                UP = Ucm + (postCollisionRelU*mQ/(mP + mQ)); // P is the NON-DISSOCIATING molecule.
                UQ = Ucm - (postCollisionRelU*mP/(mP + mQ)); // UQ is used as Ucm for atomic split.

                // P remains NON-DISSOCIATED (but internal and translational energy modified).

                label p1 = 0;
                label p2 = 1;

                const label& typeId1 = productsToDiss_[1][p1]; //note [1] corresponds to the second set of products in chemReactDict
                const label& typeId2 = productsToDiss_[1][p2];

                //Mass of Product one and two
                scalar mP1 = cloud_.constProps(typeId1).mass();
                scalar mP2 = cloud_.constProps(typeId2).mass();

                // P remains NON-DISSOCIATED (but internal and translational energy modified).
                p.vibLevel() = EVibP/(cloud_.constProps(typeIdP).thetaV()*physicoChemical::k.value());
                p.ERot() = ERotP;
                p.U() = UP;  
               
                scalar mRatoms = mP1*mP2/(mP1 + mP2);

                vector UcmAtomsQ = UQ;

//                 scalar cRatoms = sqrt(2.0*Ecplx/mRatoms);
                
                scalar cRatoms = sqrt(2.0*ERotQ/mRatoms);

                // Variable Hard Sphere collision part for atoms
            
                cosTheta = 2.0*cloud_.rndGen().sample01<scalar>() - 1.0;
            
                sinTheta = sqrt(1.0 - cosTheta*cosTheta);
            
                phi = twoPi*cloud_.rndGen().sample01<scalar>();
            
                vector postCollisionRelU2 = cRatoms
                *vector
                    (
                        cosTheta,
                        sinTheta*cos(phi),
                        sinTheta*sin(phi)
                    );

                vector uP1 = UcmAtomsQ + (postCollisionRelU2*mP2/(mP1 + mP2));
                vector uP2 = UcmAtomsQ - (postCollisionRelU2*mP1/(mP1 + mP2));

                // Molecule Q will dissociate into 2 atoms.
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
                q.vibLevel() = 0;
                q.ERot() = 0.0;
                
                label classificationQ = q.classification();
                scalar RWF = q.RWF();
                
                // insert new product 2
                cloud_.addNewParcel
                (
                    position,
                    uP2,
                    RWF,
                    0.0,
                    0,
                    0,
                    cell,
                    tetFace,
                    tetPt,
                    typeId2,
                    -1,
                    classificationQ
                );
            }
        }

        if(dissoP == true && dissoQ == true)
        {
            nABReactions_++;
            nCDReactions_++;

            relax_ = false;
    
            if(allowSplitting_)
            {
                scalar heatOfReactionJoules = (heatOfReactionAB_+heatOfReactionCD_)*physicoChemical::k.value();
                
                //Both dissociating, so no redistribution to vibrational modes
                
                translationalEnergy = translationalEnergy + heatOfReactionJoules + EVibP + EVibQ;
                
//                 if (0.2 > cloud_.rndGen().sample01<scalar>())
//                 {
//                     scalar EcP = translationalEnergy + ERotP;
//                     
//                     scalar energyRatio = 0.0;
//                     
//                     if(rotationalDofP == 2.0)
//                     {
//                         energyRatio = 1.0 - pow(cloud_.rndGen().sample01<scalar>(),(1.0/ChiB));
//                     }
//                     else
//                     {
//                         scalar ChiA = 0.5*rotationalDofP;
//                         
//                         energyRatio = cloud_.energyRatio(ChiA, ChiB);
//                     }
// 
//                     ERotP = energyRatio*EcP;
//                 
//                     translationalEnergy = EcP - ERotP;
//                 }
//                 
//                 //Molecule Q is dissociating, so no vibrational relaxation to it
//                 
//                 if (0.2 > cloud_.rndGen().sample01<scalar>())
//                 {
//                     scalar EcQ = translationalEnergy + ERotQ;
//                     
//                     scalar energyRatio = 0.0;
//                     
//                     if(rotationalDofQ == 2.0)
//                     {
//                         energyRatio = 1.0 - pow(cloud_.rndGen().sample01<scalar>(),(1.0/ChiB));
//                     }
//                     else
//                     {
//                         scalar ChiA = 0.5*rotationalDofQ;
//                         
//                         energyRatio = cloud_.energyRatio(ChiA, ChiB);
//                     }
// 
//                     ERotQ = energyRatio*EcQ;
//                 
//                     translationalEnergy = EcQ - ERotQ;
//                 }
                
                scalar relVelNonDissoMol = sqrt(2.0*translationalEnergy/mR);

                //center of mass velocity of all reactant particles
                vector Ucm = (mP*UP + mQ*UQ)/(mP + mQ);

                // Variable Hard Sphere collision part for collision of molecules
        
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

                label p1 = 0;
                label p2 = 1;

                const label& typeId1 = productsToDiss_[0][p1]; //note [0] corresponds to the 1st set of products in chemReactDict
                const label& typeId2 = productsToDiss_[0][p2];

                const label& typeId3 = productsToDiss_[1][p1]; //note [1] corresponds to the 2nd set of products in chemReactDict
                const label& typeId4 = productsToDiss_[1][p2];

                //Mass of Product one and two
                scalar mP1 = cloud_.constProps(typeId1).mass();
                scalar mP2 = cloud_.constProps(typeId2).mass();
                scalar mP3 = cloud_.constProps(typeId3).mass();
                scalar mP4 = cloud_.constProps(typeId4).mass();


                scalar mRatoms1 = mP1*mP2/(mP1 + mP2);
                scalar mRatoms2 = mP3*mP4/(mP3 + mP4);
                
                scalar cRatoms1 = sqrt(2.0*ERotP/mRatoms1);
                scalar cRatoms2 = sqrt(2.0*ERotQ/mRatoms2);

                // Variable Hard Sphere collision part for atoms1
            
                cosTheta = 2.0*cloud_.rndGen().sample01<scalar>() - 1.0;
            
                sinTheta = sqrt(1.0 - cosTheta*cosTheta);
            
                phi = twoPi*cloud_.rndGen().sample01<scalar>();
            
                vector postCollisionRelU1 = cRatoms1
                *vector
                    (
                        cosTheta,
                        sinTheta*cos(phi),
                        sinTheta*sin(phi)
                    );

                vector uP1 = UP + (postCollisionRelU1*mP2/(mP1 + mP2));
                vector uP2 = UP - (postCollisionRelU1*mP1/(mP1 + mP2));

                // Variable Hard Sphere collision part for atoms2
            
                cosTheta = 2.0*cloud_.rndGen().sample01<scalar>() - 1.0;
            
                sinTheta = sqrt(1.0 - cosTheta*cosTheta);
            
                phi = twoPi*cloud_.rndGen().sample01<scalar>();
            
                vector postCollisionRelU2 = cRatoms2
                *vector
                    (
                        cosTheta,
                        sinTheta*cos(phi),
                        sinTheta*sin(phi)
                    );

                vector uQ1 = UQ + (postCollisionRelU2*mP4/(mP3 + mP4));
                vector uQ2 = UQ - (postCollisionRelU2*mP3/(mP3 + mP4));

                vector position = p.position();
                
                label cellP = -1;
                label tetFaceP = -1;
                label tetPtP = -1;

                mesh_.findCellFacePt
                (
                    position,
                    cellP,
                    tetFaceP,
                    tetPtP
                );
                
                p.typeId() = typeId1;
                p.U() = uP1;
                p.vibLevel() = 0;
                p.ERot() = 0.0;
                
                label classificationP = p.classification();
                scalar RWF = p.RWF();

                // insert new product 2
                cloud_.addNewParcel
                (
                    position,
                    uP2,
                    RWF,
                    0.0,
                    0,
                    0,
                    cellP,
                    tetFaceP,
                    tetPtP,
                    typeId2,
                    -1,
                    classificationP
                );

                // Molecule Q will dissociate into 2 atoms.
                
                position = q.position();
                
                label cellQ = -1;
                label tetFaceQ = -1;
                label tetPtQ = -1;

                mesh_.findCellFacePt
                (
                    position,
                    cellQ,
                    tetFaceQ,
                    tetPtQ
                );
                
                q.typeId() = typeId3;
                q.U() = uQ1;
                q.vibLevel() = 0;
                q.ERot() = 0.0;
                
                label classificationQ = q.classification();
                scalar RWF2 = q.RWF();

                // insert new product 4
                cloud_.addNewParcel
                (
                    position,
                    uQ2,
                    RWF2,
                    0.0,
                    0,
                    0,
                    cellQ,
                    tetFaceQ,
                    tetPtQ,
                    typeId4,
                    -1,
                    classificationQ
                );
            }
        }
        
        if(dissoP == false && dissoQ == false)
        {
            relax_ = true;
        }
    }

//////////////////////////////////////////////////////////////////////////////////////////////
  
    if(typeIdP == reactantIds_[1] && typeIdQ == reactantIds_[0])
    {
        vector UP = p.U();
        vector UQ = q.U();
        scalar ERotP = p.ERot();
        scalar ERotQ = q.ERot();
        scalar EVibP = p.vibLevel()*cloud_.constProps(typeIdP).thetaV()*physicoChemical::k.value();
        scalar EVibQ = q.vibLevel()*cloud_.constProps(typeIdQ).thetaV()*physicoChemical::k.value();

        scalar mP = cloud_.constProps(typeIdP).mass();
        scalar mQ = cloud_.constProps(typeIdQ).mass();

        scalar mR = mP*mQ/(mP + mQ);
        scalar cRsqr = magSqr(UP - UQ);
        scalar translationalEnergy = 0.5*mR*cRsqr;
        
        scalar thetaVP = cloud_.constProps(typeIdP).thetaV();
        scalar thetaVQ = cloud_.constProps(typeIdQ).thetaV();
        
        scalar thetaDP = cloud_.constProps(typeIdP).thetaD();
        scalar thetaDQ = cloud_.constProps(typeIdQ).thetaD();
        
        scalar ZrefP = cloud_.constProps(typeIdP).Zref();
        scalar ZrefQ = cloud_.constProps(typeIdQ).Zref();
        
        scalar refTempZvP = cloud_.constProps(typeIdP).TrefZv();
        scalar refTempZvQ = cloud_.constProps(typeIdQ).TrefZv();
        
        scalar rotationalDofP = cloud_.constProps(typeIdP).rotationalDegreesOfFreedom();
        scalar rotationalDofQ = cloud_.constProps(typeIdQ).rotationalDegreesOfFreedom();

        scalar omegaPQ =
        0.5
        *(
            cloud_.constProps(typeIdP).omega()
            + cloud_.constProps(typeIdQ).omega()
        );
        
        scalar ChiB = 2.5 - omegaPQ;

        bool dissoP = false;
        bool dissoQ = false;

        scalar EcPQ = 0.0;
        label idP = cloud_.constProps(typeIdP).charDissQuantumLevel();
        label imaxP = 0;

        // calculate if a dissociation of species P is possible
        EcPQ = translationalEnergy + EVibP;

        imaxP = EcPQ/(physicoChemical::k.value()*cloud_.constProps(typeIdP).thetaV());
        
        if(imaxP-idP > 0)
        {
            dissoP = true;
        }

        scalar EcQP = 0.0;
        label idQ = cloud_.constProps(typeIdQ).charDissQuantumLevel();
        label imaxQ = 0;

        // calculate if a dissociation of species Q is possible
        EcQP = translationalEnergy + EVibQ;
        
        if(dissoP == true)
        {
            EcQP += (heatOfReactionCD_*physicoChemical::k.value());
        }

        imaxQ = EcQP/(physicoChemical::k.value()*cloud_.constProps(typeIdQ).thetaV());
        
        if(imaxQ-idQ > 0) 
        {
            dissoQ = true;
        }
        
        relax_ = true;

        if(dissoP == true && dissoQ == false)        
        {
            nCDReactions_++;

            relax_ = false;
            
            if(allowSplitting_)
            {
                scalar heatOfReactionJoules = heatOfReactionCD_*physicoChemical::k.value();
                
                translationalEnergy = translationalEnergy + heatOfReactionJoules + EVibP;
                
                if (0.2 > cloud_.rndGen().sample01<scalar>())
                {
                    scalar EcP = translationalEnergy + ERotP;
                    
                    scalar energyRatio = 0.0;
                    
                    if(rotationalDofP == 2.0)
                    {
                        energyRatio = 1.0 - pow(cloud_.rndGen().sample01<scalar>(),(1.0/ChiB));
                    }
                    else
                    {
                        scalar ChiA = 0.5*rotationalDofP;
                        
                        energyRatio = cloud_.energyRatio(ChiA, ChiB);
                    }

                    ERotP = energyRatio*EcP;
                
                    translationalEnergy = EcP - ERotP;
                }
                
                scalar EcQ = translationalEnergy + EVibQ;
                label iMaxQ = (EcQ /(physicoChemical::k.value()*thetaVQ));

                if(iMaxQ > SMALL)
                {
                    // - Bird equations 5.42 and 11.34 gave this denominator
                    scalar TCollQ = (iMaxQ*thetaVQ) / (3.5 - omegaPQ); 
                    
                    scalar pow1 = pow((thetaDQ/TCollQ),0.333) - 1.0;

                    scalar pow2 = pow ((thetaDQ/refTempZvQ),0.333) -1.0;
                    
                    // - vibrational collision number (equation 2, Bird 2010)
                    scalar ZvQ1 = pow((thetaDQ/TCollQ),omegaPQ); 
                    
                    scalar ZvQ2 = pow(ZrefQ*(pow((thetaDQ/refTempZvQ),(-1.0*omegaPQ))),(pow1/pow2));
                    
                    scalar ZvQ = ZvQ1*ZvQ2;
//                     scalar ZvQ = 50.0;
                        
                    scalar inverseVibrationalCollisionNumberQ = 1.0/(5.0*ZvQ);
                
                    if(inverseVibrationalCollisionNumberQ > cloud_.rndGen().sample01<scalar>())
                    {
                        label iDashQ = 0; // post-collision quantum number
                        scalar func = 0.0;

                        do // acceptance - rejection 
                        {
                            iDashQ = cloud_.rndGen().integer(0,iMaxQ);
                            EVibQ = iDashQ*physicoChemical::k.value()*thetaVQ;
                            func = pow((1.0 - (EVibQ / EcQ)),(1.5 - omegaPQ));
                    
                        } while( !(func > cloud_.rndGen().sample01<scalar>()) );
                
                        translationalEnergy = EcQ - EVibQ;
                    }
                }
                
                if (0.2 > cloud_.rndGen().sample01<scalar>())
                {
                    scalar EcQ = translationalEnergy + ERotQ;
                    
                    scalar energyRatio = 0.0;
                    
                    if(rotationalDofQ == 2.0)
                    {
                        energyRatio = 1.0 - pow(cloud_.rndGen().sample01<scalar>(),(1.0/ChiB));
                    }
                    else
                    {
                        scalar ChiA = 0.5*rotationalDofQ;
                        
                        energyRatio = cloud_.energyRatio(ChiA, ChiB);
                    }

                    ERotQ = energyRatio*EcQ;
                
                    translationalEnergy = EcQ - ERotQ;
                }
                
                scalar relVelNonDissoMol = sqrt((2.0*translationalEnergy)/mR);

                //center of mass velocity of all reactant particles
                vector Ucm = (mP*UP + mQ*UQ)/(mP + mQ);

                // Variable Hard Sphere collision part for collision of molecules
        
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
    
                UP = Ucm + (postCollisionRelU*mQ/(mP + mQ)); // UP is used as Ucm for atomic split.
                UQ = Ucm - (postCollisionRelU*mP/(mP + mQ)); // Q is the NON-DISSOCIATING molecule.

                // Q remains NON-DISSOCIATED (but internal and translational energy modified).

                label p1 = 0;
                label p2 = 1;

                const label& typeId1 = productsToDiss_[1][p1]; //note [1] corresponds to the second set of products in chemReactDict
                const label& typeId2 = productsToDiss_[1][p2];

                //Mass of Product one and two
                scalar mP1 = cloud_.constProps(typeId1).mass();
                scalar mP2 = cloud_.constProps(typeId2).mass();

                // Q remains NON-DISSOCIATED (but internal and translational energy modified).
                q.vibLevel() = EVibQ/(cloud_.constProps(typeIdQ).thetaV()*physicoChemical::k.value());
                q.ERot() = ERotQ;
                q.U() = UQ;

                scalar mRatoms = mP1*mP2/(mP1 + mP2);

                vector UcmAtomsP = UP;
                
                scalar cRatoms = sqrt(2.0*ERotP/mRatoms);

                // Variable Hard Sphere collision part for atoms
            
                cosTheta = 2.0*cloud_.rndGen().sample01<scalar>() - 1.0;
            
                sinTheta = sqrt(1.0 - cosTheta*cosTheta);
            
                phi = twoPi*cloud_.rndGen().sample01<scalar>();
            
                vector postCollisionRelU2 = cRatoms
                *vector
                    (
                        cosTheta,
                        sinTheta*cos(phi),
                        sinTheta*sin(phi)
                    );

                vector uP1 = UcmAtomsP + (postCollisionRelU2*mP2/(mP1 + mP2));
                vector uP2 = UcmAtomsP - (postCollisionRelU2*mP1/(mP1 + mP2));

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
                p.vibLevel() = 0;
                p.ERot() = 0.0;
                
                label classificationP = p.classification();
                scalar RWF = p.RWF();
                
                // insert new product 2
                cloud_.addNewParcel
                (
                    position,
                    uP2,
                    RWF,
                    0.0,
                    0,
                    0,
                    cell,
                    tetFace,
                    tetPt,
                    typeId2,
                    -1,
                    classificationP
                );
            }
        }

        if(dissoQ == true && dissoP == false)
        {
            nABReactions_++;
            relax_ = false; 
            
            if(allowSplitting_)
            {
                scalar heatOfReactionJoules = heatOfReactionAB_*physicoChemical::k.value();
                                
                translationalEnergy = translationalEnergy + heatOfReactionJoules + EVibQ;   
                
//                 scalar EcP = translationalEnergy + EVibP;
//                 label iMaxP = (EcP /(physicoChemical::k.value()*thetaVP));
// 
//                 if(iMaxP > SMALL)
//                 {
//                     // - Bird equations 5.42 and 11.34 gave this denominator
//                     scalar TCollP = (iMaxP*thetaVP) / (3.5 - omegaPQ); 
//                     
//                     scalar pow1 = pow((thetaDP/TCollP),0.333) - 1.0;
// 
//                     scalar pow2 = pow ((thetaDP/refTempZvP),0.333) -1.0;
//                     
//                     // - vibrational collision number (equation 2, Bird 2010)
//                     scalar ZvP1 = pow((thetaDP/TCollP),omegaPQ); 
//                     
//                     scalar ZvP2 = pow(ZrefP*(pow((thetaDP/refTempZvP),(-1.0*omegaPQ))),(pow1/pow2));
//                     
//                     scalar ZvP = ZvP1*ZvP2;
// //                     scalar ZvP = 50.0;
//                         
//                     scalar inverseVibrationalCollisionNumberQ = 1.0/ZvP;
//                 
//                     if(inverseVibrationalCollisionNumberQ > cloud_.rndGen().sample01<scalar>())
//                     {
//                         label iDashP = 0; // post-collision quantum number
//                         scalar func = 0.0;
// 
//                         do // acceptance - rejection 
//                         {
//                             iDashP = cloud_.rndGen().integer(0,iMaxP);
//                             EVibP = iDashP*physicoChemical::k.value()*thetaVP;
//                             func = pow((1.0 - (EVibP / EcP)),(1.5 - omegaPQ));
//                     
//                         } while( !(func > cloud_.rndGen().sample01<scalar>()) );
//                 
//                         translationalEnergy = EcP - EVibP;
//                     }
//                 }
//                 
//                 if (0.2 > cloud_.rndGen().sample01<scalar>())
//                 {
//                     scalar EcP = translationalEnergy + ERotP;
//                     
//                     scalar energyRatio = 0.0;
//                     
//                     if(rotationalDofP == 2.0)
//                     {
//                         energyRatio = 1.0 - pow(cloud_.rndGen().sample01<scalar>(),(1.0/ChiB));
//                     }
//                     else
//                     {
//                         scalar ChiA = 0.5*rotationalDofP;
//                         
//                         energyRatio = cloud_.energyRatio(ChiA, ChiB);
//                     }
// 
//                     ERotP = energyRatio*EcP;
//                 
//                     translationalEnergy = EcP - ERotP;
//                 }
                
                scalar relVelNonDissoMol = sqrt((2.0*translationalEnergy)/mR);

                //center of mass velocity of all reactant particles
                vector Ucm = (mP*UP + mQ*UQ)/(mP + mQ);

                // Variable Hard Sphere collision part for collision of molecules
        
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
        
                UP = Ucm + (postCollisionRelU*mQ/(mP + mQ)); // P is the NON-DISSOCIATING molecule.

                UQ = Ucm - (postCollisionRelU*mP/(mP + mQ)); // UQ is used as Ucm for atomic split.

                // P remains NON-DISSOCIATED (but internal and translational energy modified).

                label p1 = 0;
                label p2 = 1;

                const label& typeId1 = productsToDiss_[0][p1]; //note [0] corresponds to the 1st set of products in chemReactDict
                const label& typeId2 = productsToDiss_[0][p2];

                //Mass of Product one and two
                scalar mP1 = cloud_.constProps(typeId1).mass();
                scalar mP2 = cloud_.constProps(typeId2).mass();

                // P remains NON-DISSOCIATED (but internal and translational energy modified).

                p.vibLevel() = EVibP/(cloud_.constProps(typeIdP).thetaV()*physicoChemical::k.value());
                p.ERot() = ERotP;
                p.U() = UP;

                scalar mRatoms = mP1*mP2/(mP1 + mP2);

                vector UcmAtomsQ = UQ;
                
                scalar cRatoms = sqrt(2.0*ERotQ/mRatoms);

                // Variable Hard Sphere collision part for atoms
            
                cosTheta = 2.0*cloud_.rndGen().sample01<scalar>() - 1.0;
            
                sinTheta = sqrt(1.0 - cosTheta*cosTheta);
            
                phi = twoPi*cloud_.rndGen().sample01<scalar>();
            
                vector postCollisionRelU2 = cRatoms
                *vector
                    (
                        cosTheta,
                        sinTheta*cos(phi),
                        sinTheta*sin(phi)
                    );

                vector uP1 = UcmAtomsQ + (postCollisionRelU2*mP2/(mP1 + mP2));
                vector uP2 = UcmAtomsQ - (postCollisionRelU2*mP1/(mP1 + mP2));

                // Molecule Q will dissociate into 2 atoms.
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
                q.vibLevel() = 0;
                q.ERot() = 0.0;
                
                label classificationQ = q.classification();
                scalar RWF = q.RWF();
                
                // insert new product 2
                cloud_.addNewParcel
                (
                    position,
                    uP2,
                    RWF,
                    0.0,
                    0,
                    0,
                    cell,
                    tetFace,
                    tetPt,
                    typeId2,
                    -1,
                    classificationQ
                );
            }
        } 

        if(dissoP == true && dissoQ == true)
        {
            nABReactions_++;
            nCDReactions_++;

            relax_ = false;
    
            if(allowSplitting_)
            {
                scalar heatOfReactionJoules = (heatOfReactionCD_+heatOfReactionAB_)*physicoChemical::k.value();
                
                translationalEnergy = translationalEnergy + heatOfReactionJoules + EVibP + EVibQ;
                
                scalar relVelNonDissoMol = sqrt(2.0*translationalEnergy/mR);

                //center of mass velocity of all reactant particles
                vector Ucm = (mP*UP + mQ*UQ)/(mP + mQ);

                // Variable Hard Sphere collision part for collision of molecules
        
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

                label p1 = 0;
                label p2 = 1;

                const label& typeId1 = productsToDiss_[1][p1]; //note [0] corresponds to the 1st set of products in chemReactDict
                const label& typeId2 = productsToDiss_[1][p2];

                const label& typeId3 = productsToDiss_[0][p1]; //note [1] corresponds to the 2nd set of products in chemReactDict
                const label& typeId4 = productsToDiss_[0][p2];

                //Mass of Product one and two
                scalar mP1 = cloud_.constProps(typeId1).mass();
                scalar mP2 = cloud_.constProps(typeId2).mass();
                scalar mP3 = cloud_.constProps(typeId3).mass();
                scalar mP4 = cloud_.constProps(typeId4).mass();


                scalar mRatoms1 = mP1*mP2/(mP1 + mP2);
                scalar mRatoms2 = mP3*mP4/(mP3 + mP4);
                
                scalar cRatoms1 = sqrt(2.0*ERotP/mRatoms1);
                scalar cRatoms2 = sqrt(2.0*ERotQ/mRatoms2);

                // Variable Hard Sphere collision part for atoms1
            
                cosTheta = 2.0*cloud_.rndGen().sample01<scalar>() - 1.0;
            
                sinTheta = sqrt(1.0 - cosTheta*cosTheta);
            
                phi = twoPi*cloud_.rndGen().sample01<scalar>();
            
                vector postCollisionRelU1 = cRatoms1
                *vector
                    (
                        cosTheta,
                        sinTheta*cos(phi),
                        sinTheta*sin(phi)
                    );

                vector uP1 = UP + (postCollisionRelU1*mP2/(mP1 + mP2));
                vector uP2 = UP - (postCollisionRelU1*mP1/(mP1 + mP2));

                // Variable Hard Sphere collision part for atoms2
            
                cosTheta = 2.0*cloud_.rndGen().sample01<scalar>() - 1.0;
            
                sinTheta = sqrt(1.0 - cosTheta*cosTheta);
            
                phi = twoPi*cloud_.rndGen().sample01<scalar>();
            
                vector postCollisionRelU2 = cRatoms2
                *vector
                    (
                        cosTheta,
                        sinTheta*cos(phi),
                        sinTheta*sin(phi)
                    );

                vector uQ1 = UQ + (postCollisionRelU2*mP4/(mP3 + mP4));
                vector uQ2 = UQ - (postCollisionRelU2*mP3/(mP3 + mP4));


                // Molecule P will dissociate into 2 atoms.
                vector position = p.position();
                
                label cellP = -1;
                label tetFaceP = -1;
                label tetPtP = -1;

                mesh_.findCellFacePt
                (
                    position,
                    cellP,
                    tetFaceP,
                    tetPtP
                );
                
                p.typeId() = typeId1;
                p.U() = uP1;
                p.vibLevel() = 0;
                p.ERot() = 0.0;
                
                label classificationP = p.classification();
                scalar RWF = p.RWF();

                // insert new product 2
                cloud_.addNewParcel
                (
                    position,
                    uP2,
                    RWF,
                    0.0,
                    0,
                    0,
                    cellP,
                    tetFaceP,
                    tetPtP,
                    typeId2,
                    -1,
                    classificationP
                );

                // Molecule Q will dissociate into 2 atoms.
                position = q.position();
                
                label cellQ = -1;
                label tetFaceQ = -1;
                label tetPtQ = -1;

                mesh_.findCellFacePt
                (
                    position,
                    cellQ,
                    tetFaceQ,
                    tetPtQ
                );
                
                q.typeId() = typeId3;
                q.U() = uQ1;
                q.vibLevel() = 0;
                q.ERot() = 0.0;
                
                label classificationQ = q.classification();
                scalar RWF2 = q.RWF();

                // insert new product 4
                cloud_.addNewParcel
                (
                    position,
                    uQ2,
                    RWF2,
                    0.0,
                    0,
                    0,
                    cellQ,
                    tetFaceQ,
                    tetPtQ,
                    typeId4,
                    -1,
                    classificationQ
                );
            }
        }
        
        if(dissoQ == false && dissoP == false)
        {
            relax_ = true;
        }
    }   
}
//&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
void  dissociationTypeIDissimilarSpecies::outputResults(const label& counterIndex)
{  
    if(writeRatesToTerminal_ == true)
    {
            // measure density 
//         if(counterIndex == 1)
//         {
            const List< DynamicList<dsmcParcel*> >& cellOccupancy
                = cloud_.cellOccupancy();

            List<label> mols;
	mols.append(0); mols.append(0);

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
            label nABReactions = nABReactions_;
            label nCDReactions = nCDReactions_;

            //- Parallel communication
            if(Pstream::parRun())
            {
                reduce(volume, sumOp<scalar>());
                reduce(mols[0], sumOp<label>());
                reduce(mols[1], sumOp<label>());
                reduce(nABReactions, sumOp<label>());
                reduce(nCDReactions, sumOp<label>());
            }

            numberDensities_[0] = (mols[0]*cloud().nParticle())/volume;
            numberDensities_[1] = (mols[1]*cloud().nParticle())/volume;
//         }

        const scalar& deltaT = mesh_.time().deltaT().value();

        word reactantMolA = cloud_.typeIdList()[reactantIds_[0]];
        word reactantMolB = cloud_.typeIdList()[reactantIds_[1]];
        
        word dissociationProductMolA = cloud_.typeIdList()[productsToDiss_[0][0]];
        word dissociationProductMolB = cloud_.typeIdList()[productsToDiss_[0][1]];
        word dissociationProductMolC = cloud_.typeIdList()[productsToDiss_[1][0]];
        word dissociationProductMolD = cloud_.typeIdList()[productsToDiss_[1][1]];

        if((numberDensities_[0] > 0.0) && (numberDensities_[1] > 0.0))
        {
            scalar reactionRate1 = 0.0;
            scalar reactionRate2 = 0.0;

            reactionRate1 =
            (
                nABReactions
                * cloud_.nParticle()
            )/(counterIndex*deltaT*numberDensities_[0]* numberDensities_[1]*volume);
        
            reactionRate2 =
            (
                nCDReactions
                * cloud_.nParticle()
            )/(counterIndex*deltaT*numberDensities_[0]* numberDensities_[1]*volume);

            Info << "Dissociation type I reaction " <<  reactantMolA << " + " << reactantMolB << " --> " <<
            dissociationProductMolA << " + " << dissociationProductMolB << " + " << reactantMolB <<
            ", reaction rate = " << reactionRate1  << nl      
            << "Dissociation type I reaction " <<  reactantMolB << " + " << reactantMolA << " --> " <<
            dissociationProductMolC << " + " << dissociationProductMolD << " + " << reactantMolA <<
            ", reaction rate = " << reactionRate2
            << endl;
        }
    }
    else
    {
        scalar nABReactions = nABReactions_;
        scalar nCDReactions = nCDReactions_;
        
        //- Parallel communication
        if(Pstream::parRun())
        {
            reduce(nABReactions, sumOp<scalar>());
            reduce(nCDReactions, sumOp<scalar>());
        }
        
        if(nABReactions > VSMALL)
        {
            word reactantMolA = cloud_.typeIdList()[reactantIds_[0]];
            word reactantMolB = cloud_.typeIdList()[reactantIds_[1]];
            
            word dissociationProductMolA = cloud_.typeIdList()[productsToDiss_[0][0]];
            word dissociationProductMolB = cloud_.typeIdList()[productsToDiss_[0][1]];
            
            Info << "Dissociation type I reaction " <<  reactantMolA << " + " << reactantMolB << " --> " <<
                dissociationProductMolA << " + " << dissociationProductMolB << " + " << reactantMolB <<
                " is active." << endl;
        }
        
        if(nCDReactions > VSMALL)
        {
            word reactantMolA = cloud_.typeIdList()[reactantIds_[0]];
            word reactantMolB = cloud_.typeIdList()[reactantIds_[1]];
            
            word dissociationProductMolC = cloud_.typeIdList()[productsToDiss_[1][0]];
            word dissociationProductMolD = cloud_.typeIdList()[productsToDiss_[1][1]];
               
            Info << "Dissociation type I reaction " <<  reactantMolB << " + " << reactantMolA << " --> " <<
                dissociationProductMolC << " + " << dissociationProductMolD << " + " << reactantMolA <<
                " is active." << endl;
        }

    }

    nReactionsPerTimeStep_ = 0.0;

}


const bool& dissociationTypeIDissimilarSpecies::relax() const
{
    return relax_;
}

}
// End namespace Foam

// ************************************************************************* //
