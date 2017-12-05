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

#include "mixedTypeIIDissociationReverseExchange.H"
#include "addToRunTimeSelectionTable.H"
#include "fvc.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(mixedTypeIIDissociationReverseExchange, 0);

addToRunTimeSelectionTable(dsmcReaction, mixedTypeIIDissociationReverseExchange, dictionary);



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
mixedTypeIIDissociationReverseExchange::mixedTypeIIDissociationReverseExchange
(
    Time& t,
    dsmcCloud& cloud,
    const dictionary& dict
)
:
    dsmcReaction(t, cloud, dict),
    propsDict_(dict.subDict(typeName + "Properties")),
    reactantIds_(),
    exchangeProductIds_(),
    dissociationProductIds_(),
    activationEnergy_(readScalar(propsDict_.lookup("activationEnergy"))),
    heatOfReactionDiss_(readScalar(propsDict_.lookup("heatOfReactionDiss"))),
    heatOfReactionExch_(readScalar(propsDict_.lookup("heatOfReactionExch"))),
    aCoeff_(readScalar(propsDict_.lookup("aCoeff"))),
    bCoeff_(readScalar(propsDict_.lookup("bCoeff"))),
    nTotExchangeReactions_(0),
    reactionName_(propsDict_.lookup("reactionName")),
    relax_(true),
    allowSplitting_(false),
    writeRatesToTerminal_(false),
    volume_(0.0),
    numberDensities_(2, 0.0)
{

}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

mixedTypeIIDissociationReverseExchange::~mixedTypeIIDissociationReverseExchange()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void mixedTypeIIDissociationReverseExchange::initialConfiguration()
{
    setProperties();
}

void mixedTypeIIDissociationReverseExchange::setProperties()
{
    // reading in reactants
    
//     reactionName_(propsDict_.lookup("reactionName"));

    const List<word> reactantMolecules (propsDict_.lookup("reactantMolecules"));

    if(reactantMolecules.size() != 2)
    {
        FatalErrorIn("mixedTypeIIDissociationReverseExchange::setProperties()")
            << "There should be two or more reactants, instead of " 
            << reactantMolecules.size() << nl 
            << exit(FatalError);
    }
    
    if(reactantMolecules[0] == reactantMolecules[1])
    {
        FatalErrorIn("mixedTypeIIDissociationReverseExchange::setProperties()")
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
            FatalErrorIn("mixedTypeIIDissociationReverseExchange::setProperties()")
                << "Cannot find type id: " << reactantMolecules[r] << nl 
                << exit(FatalError);
        }
    }
    
    // check that the first reactants is a 'MOLECULE' 

    const scalar& rDof1 = cloud_.constProps(reactantIds_[0]).rotationalDegreesOfFreedom();

    if(rDof1 < 1)
    {
        FatalErrorIn("mixedTypeIIDissociationReverseExchange::setProperties()")
            << "Reactant 1 must be a molecule (not an atom): " << reactantMolecules[0] 
            << nl 
            << exit(FatalError);
    }
    
    // check that the second reactant is an 'ATOM' 

    const scalar& rDof2 = cloud_.constProps(reactantIds_[1]).rotationalDegreesOfFreedom();

    if(rDof2 > 0)
    {
        FatalErrorIn("mixedTypeIIDissociationReverseExchange::setProperties()")
            << "Reactant 2 must be an atom (not a molecule): " << reactantMolecules[1] 
            << nl 
            << exit(FatalError);
    }

    // reading in products

    const List<word> exchangeProductMolecules (propsDict_.lookup("productsOfExchangeReaction"));

    if(exchangeProductMolecules.size() != 2)
    {
        FatalErrorIn("mixedTypeIIDissociationReverseExchange::setProperties()")
            << "There should be two mixedTypeIIDissociationReverseExchange reaction products, instead of " 
            << reactantMolecules.size() << nl 
            << exit(FatalError);
    }
    
    if(exchangeProductMolecules[0] == exchangeProductMolecules[1])
    {
        FatalErrorIn("mixedTypeIIDissociationReverseExchange::setProperties()")
            << "Exchange reaction product molecules cannot be same species." << nl
            << exit(FatalError);
    }

    exchangeProductIds_.setSize(exchangeProductMolecules.size(), -1);

    forAll(exchangeProductIds_, r)
    {
        exchangeProductIds_[r] = findIndex(cloud_.typeIdList(), exchangeProductMolecules[r]);

        // check that reactants belong to the typeIdList (constant/dsmcProperties)
        if(exchangeProductIds_[r] == -1)
        {
            FatalErrorIn("mixedTypeIIDissociationReverseExchange::setProperties()")
                << "Cannot find type id: " << exchangeProductMolecules[r] << nl 
                << exit(FatalError);
        }
    }
    
    const List<word> dissociationProductMolecules (propsDict_.lookup("productsOfDissociatedMolecule"));
    
    if(dissociationProductMolecules.size() != 2)
    {
        FatalErrorIn("mixedTypeIIDissociationReverseExchange::setProperties()")
            << "There should be two mixedTypeIIDissociationReverseExchange reaction products, instead of " 
            << reactantMolecules.size() << nl 
            << exit(FatalError);
    }

    dissociationProductIds_.setSize(dissociationProductMolecules.size(), -1);

    forAll(dissociationProductIds_, r)
    {
        dissociationProductIds_[r] = findIndex(cloud_.typeIdList(), dissociationProductMolecules[r]);

        // check that reactants belong to the typeIdList (constant/dsmcProperties)
        if(dissociationProductIds_[r] == -1)
        {
            FatalErrorIn("mixedTypeIIDissociationReverseExchange::setProperties()")
                << "Cannot find type id: " << dissociationProductMolecules[r] << nl 
                << exit(FatalError);
        }
    }
    
    //check that first exchange product is a 'MOLECULE' (not an 'ATOM')
    
    const scalar& rDof3 = cloud_.constProps(exchangeProductIds_[0]).rotationalDegreesOfFreedom();

    if(rDof3 < 0)
    {
        FatalErrorIn("mixedTypeIIDissociationReverseExchange::setProperties()")
            << "First product of the exchange reaction must be a molecule (not an atom): " << exchangeProductMolecules[0] 
            << nl 
            << exit(FatalError);
    }
    
    //check that second exchange product is an 'ATOM' (not a 'MOLECULE')
    
    const scalar& rDof4 = cloud_.constProps(exchangeProductIds_[1]).rotationalDegreesOfFreedom();

    if(rDof4 > 0)
    {
        FatalErrorIn("mixedTypeIIDissociationReverseExchange::setProperties()")
            << "Second product of the exchange reaction must be an atom (not a molecule): " << exchangeProductMolecules[1] 
            << nl 
            << exit(FatalError);
    }
    
        const scalar& rDof5= cloud_.constProps(dissociationProductIds_[0]).rotationalDegreesOfFreedom();

    if(rDof5 > 0)
    {
        FatalErrorIn("mixedTypeIIDissociationReverseExchange::setProperties()")
            << "First product of the dissociation reaction must be an atom (not a molecule): " << dissociationProductMolecules[0] 
            << nl 
            << exit(FatalError);
    }
    
    //check that second exchange product is an 'ATOM' (not a 'MOLECULE')
    
    const scalar& rDof6 = cloud_.constProps(dissociationProductIds_[1]).rotationalDegreesOfFreedom();

    if(rDof6 > 0)
    {
        FatalErrorIn("mixedTypeIIDissociationReverseExchange::setProperties()")
            << "Second product of the exchange reaction must be an atom (not a molecule): " << dissociationProductMolecules[1] 
            << nl 
            << exit(FatalError);
    }
    
    activationEnergy_ *= physicoChemical::k.value();
}

bool mixedTypeIIDissociationReverseExchange::tryReactMolecules(const label& typeIdP, const label& typeIdQ)
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

void mixedTypeIIDissociationReverseExchange::reaction
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


void mixedTypeIIDissociationReverseExchange::reaction
(
    dsmcParcel& p,
    dsmcParcel& q
)

{

    label typeIdP = p.typeId();
    label typeIdQ = q.typeId();

    //if particle p is the molecule and q is the atom...
    if(typeIdP == reactantIds_[0] && typeIdQ == reactantIds_[1]) // AB + C ---> A + B + C (diss type II) or AB + C ---> AC + B (forward exchange)
    {
        vector UP = p.U();
        vector UQ = q.U();
        scalar ERotP = p.ERot();
        scalar EVibP = p.vibLevel()*cloud_.constProps(typeIdP).thetaV()*physicoChemical::k.value();

        scalar mP = cloud_.constProps(typeIdP).mass();

        scalar mQ = cloud_.constProps(typeIdQ).mass();

        scalar mR = mP*mQ/(mP + mQ);
        scalar cRsqr = magSqr(UP - UQ);
        scalar translationalEnergy = 0.5*mR*cRsqr;
        
        scalar thetaV = cloud_.constProps(exchangeProductIds_[0]).thetaV();
        scalar thetaD = cloud_.constProps(exchangeProductIds_[0]).thetaD();
        scalar Zref = cloud_.constProps(exchangeProductIds_[0]).Zref();
        scalar refTempZv = cloud_.constProps(exchangeProductIds_[0]).TrefZv();
        
        scalar rotationalDofP = cloud_.constProps(typeIdP).rotationalDegreesOfFreedom();
        scalar rotationalDofExch = cloud_.constProps(exchangeProductIds_[0]).rotationalDegreesOfFreedom();

        scalar omegaPQ =
        0.5
        *(
            cloud_.constProps(typeIdP).omega()
            + cloud_.constProps(typeIdQ).omega()
        );
        
        scalar ChiB = 2.5 - omegaPQ;
        
        scalar omegaPQExch =
        0.5
        *(
            cloud_.constProps(exchangeProductIds_[0]).omega()
            + cloud_.constProps(exchangeProductIds_[1]).omega()
        );
        
        scalar ChiBExch = 2.5 - omegaPQExch;

        scalar EcPQ = 0.0;
        scalar TColl = 0.0;
        label idP = -1;
        label deltaDissoIP= 0;
        label imaxP = 0;
        label iaP = 0;

        // firstly calculate dissociation probability (0 or 1).
        EcPQ = translationalEnergy + EVibP;

        imaxP = EcPQ/(physicoChemical::k.value()*cloud_.constProps(typeIdP).thetaV());

        idP = cloud_.constProps(typeIdP).thetaD()/cloud_.constProps(typeIdP).thetaV();

        deltaDissoIP = imaxP-idP;
        
        relax_ = true;

        if(deltaDissoIP > 0)
        {
            // DISSOCIATION //

            nTotDissociationReactions_++;
            
            if (allowSplitting_)
            {
                relax_ = false; 
                
                scalar heatOfReactionDissJoules = heatOfReactionDiss_*physicoChemical::k.value();
                
                translationalEnergy = translationalEnergy + heatOfReactionDissJoules + EVibP;
                
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

                label p1 = 0;
                label p2 = 1;

                const label& typeId1 = dissociationProductIds_[p1];
                const label& typeId2 = dissociationProductIds_[p2];

                //Mass of Product one and two
                scalar mP1 = cloud_.constProps(typeId1).mass();
                scalar mP2 = cloud_.constProps(typeId2).mass();

                scalar mRatoms = mP1*mP2/(mP1 + mP2);

                //center of mass velocity of all particles

                vector UcmAtoms = UP;
                
                scalar cRatoms = sqrt(2.0*ERotP/mRatoms);

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
                p.ERot() = 0.0;
                p.vibLevel() = 0;
                
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
                    0,
                    classificationP
                );
            }
        }
        
        if(deltaDissoIP <= 0)
        {
            scalar P_exch = 0.0;

            // Next, calculate exchange probability.

            TColl = (translationalEnergy/(physicoChemical::k.value()))/(2.5 - omegaPQ);// Eq.(10) of Bird QK paper.
                
            scalar heatOfReactionExchJoules = heatOfReactionExch_*physicoChemical::k.value();

            scalar activationEnergy = activationEnergy_ + (aCoeff_*pow((TColl/273.0) , bCoeff_) * fabs(heatOfReactionExchJoules));

            if(EcPQ > activationEnergy) // i.e. exchange reaction can possibly occur.
            {
                scalar summation = 0.0; // declare "summation" term.

                if(activationEnergy < cloud_.constProps(typeIdP).thetaV()*physicoChemical::k.value())
                {
                    summation = 1.0; // this refers to the first sentence in Bird's QK paper after Eq.(12).
                }
                else
                {
                    iaP = EcPQ / (physicoChemical::k.value()*cloud_.constProps(typeIdP).thetaV());

                    for(label i = 0 ; i <= iaP ; i++)
                    {
                        summation += pow((1.0 - ((i*physicoChemical::k.value()*cloud_.constProps(typeIdP).thetaV())/EcPQ)),(1.5 - omegaPQ));
                    }
                }

                P_exch = pow((1.0 - (activationEnergy/EcPQ)),(1.5 - omegaPQ))/summation;// now based on modified activation energy.
            }

            if(P_exch > cloud_.rndGen().scalar01())            
            {
            // EXCHANGE REACTION //

                nTotExchangeReactions_++;
                
                if (allowSplitting_)
                {
                    relax_ = false;
                    
                    //center of mass velocity (pre-exchange)
                    vector Ucm = (mP*UP + mQ*UQ)/(mP + mQ);

                    const label& typeIdMol = exchangeProductIds_[0];
                    const label& typeIdAtom = exchangeProductIds_[1];

                    // change species properties
                    
                    scalar mPExch = cloud_.constProps(typeIdAtom).mass();
                    scalar mQExch = cloud_.constProps(typeIdMol).mass();

                    scalar mRExch = mPExch*mQExch/(mPExch + mQExch);
                    
                    translationalEnergy = translationalEnergy + ERotP + EVibP + heatOfReactionExchJoules;
                    
                    q.typeId() = typeIdMol; // q is originally the atom, becomes the molecule
                    p.typeId() = typeIdAtom; // p is the originally the molecule, becomes the atom
                                        
                    scalar EVibQ = 0.0;
                    
                    scalar EcQ = translationalEnergy + EVibQ;
                    label iMaxQ = (EcQ /(physicoChemical::k.value()*thetaV));

                    if(iMaxQ > SMALL)
                    {
                        // - Bird equations 5.42 and 11.34 gave this denominator
                        scalar TCollQ = (iMaxQ*thetaV) / (3.5 - omegaPQExch); 
                        
                        scalar pow1 = pow((thetaD/TCollQ),0.333) - 1.0;

                        scalar pow2 = pow ((thetaD/refTempZv),0.333) -1.0;
                        
                        // - vibrational collision number (equation 2, Bird 2010)
                        scalar ZvQ1 = pow((thetaD/TCollQ),omegaPQExch); 
                        
                        scalar ZvQ2 = pow(Zref*(pow((thetaD/refTempZv),(-1.0*omegaPQExch))),(pow1/pow2));
                        
                        scalar ZvQ = ZvQ1*ZvQ2;
//                         scalar ZvQ = 50.0;
                            
                        scalar inverseVibrationalCollisionNumberQ = 1.0/ZvQ;
                    
                        if(inverseVibrationalCollisionNumberQ > cloud_.rndGen().scalar01())
                        {
                            label iDashQ = 0; // post-collision quantum number
                            scalar func = 0.0;

                            do // acceptance - rejection 
                            {
                                iDashQ = cloud_.rndGen().integer(0,iMaxQ);
                                EVibQ = iDashQ*physicoChemical::k.value()*thetaV;
                                func = pow((1.0 - (EVibQ / EcQ)),(1.5 - omegaPQExch));
                        
                            } while( !(func > cloud_.rndGen().scalar01()) );
                    
                            translationalEnergy = EcQ - EVibQ;
                        }
                    }
                    
                    scalar ERotQ = 0.0;
                    
                    if (0.2 > cloud_.rndGen().scalar01())
                    {
                        scalar EcQ = translationalEnergy + ERotQ;
                        
                        scalar energyRatio = 0.0;
                        
                        if(rotationalDofExch == 2.0)
                        {
                            energyRatio = 1.0 - pow(cloud_.rndGen().scalar01(),(1.0/ChiBExch));
                        }
                        else
                        {
                            scalar ChiA = 0.5*rotationalDofExch;
                            
                            energyRatio = cloud_.energyRatio(ChiA, ChiBExch);
                        }

                        ERotQ = energyRatio*EcQ;
                    
                        translationalEnergy = EcQ - ERotQ;
                    }
                    
                    scalar relVelExchMol = sqrt((2.0*translationalEnergy)/mRExch);

                    // Variable Hard Sphere collision part for collision of molecules
            
                    scalar cosTheta = 2.0*cloud_.rndGen().scalar01() - 1.0;
                
                    scalar sinTheta = sqrt(1.0 - cosTheta*cosTheta);
                
                    scalar phi = twoPi*cloud_.rndGen().scalar01();
                
                    vector postCollisionRelU =
                        relVelExchMol
                        *vector
                            (
                                cosTheta,
                                sinTheta*cos(phi),
                                sinTheta*sin(phi)
                            );
            
                    UP = Ucm + (postCollisionRelU*mQExch/(mPExch + mQExch)); // P changes from mol to atom.
                    UQ = Ucm - (postCollisionRelU*mPExch/(mPExch + mQExch)); // Q changes from atom to mol.

                    q.vibLevel() = EVibQ/(cloud_.constProps(typeIdQ).thetaV()*physicoChemical::k.value());
                    q.ERot() = ERotQ;
                    q.U() = UQ;
                    
                    p.U() = UP;
                    p.ERot() = 0.0; // remove p's internal energies as it's now an atom
                    p.vibLevel() = 0;
                }
            }
        }
    }
     
    // if q is the molecule and p is the atom... 
    if(typeIdP == reactantIds_[1] && typeIdQ == reactantIds_[0])// C + AB ---> A + B + C (diss type II) or C + AB ---> AC + B (forward exchange)
    {
        vector UP = p.U();
        vector UQ = q.U();
        scalar ERotQ = q.ERot();
        scalar EVibQ = q.vibLevel()*cloud_.constProps(typeIdQ).thetaV()*physicoChemical::k.value();

        scalar mP = cloud_.constProps(typeIdP).mass();
        scalar mQ = cloud_.constProps(typeIdQ).mass();

        scalar mR = mP*mQ/(mP + mQ);
        scalar cRsqr = magSqr(UP - UQ);
        scalar translationalEnergy = 0.5*mR*cRsqr;
        
        scalar thetaV = cloud_.constProps(exchangeProductIds_[0]).thetaV();
        scalar thetaD = cloud_.constProps(exchangeProductIds_[0]).thetaD();
        scalar Zref = cloud_.constProps(exchangeProductIds_[0]).Zref();
        scalar refTempZv = cloud_.constProps(exchangeProductIds_[0]).TrefZv();

        scalar rotationalDofQ = cloud_.constProps(typeIdQ).rotationalDegreesOfFreedom();
        scalar rotationalDofExch = cloud_.constProps(exchangeProductIds_[0]).rotationalDegreesOfFreedom();

        scalar omegaPQ =
        0.5
        *(
            cloud_.constProps(typeIdP).omega()
            + cloud_.constProps(typeIdQ).omega()
        );
        
        scalar ChiB = 2.5 - omegaPQ;
        
        scalar omegaPQExch =
        0.5
        *(
            cloud_.constProps(exchangeProductIds_[0]).omega()
            + cloud_.constProps(exchangeProductIds_[1]).omega()
        );
        
        scalar ChiBExch = 2.5 - omegaPQExch;

        scalar EcPQ = 0.0;
        scalar TColl = 0.0;
        label idQ = -1;
        label deltaDissoIQ= 0;
        label imaxQ = 0;
        label iaQ = 0;

        // firstly calculate dissociation probability (0 or 1).
        EcPQ = translationalEnergy + EVibQ;

        imaxQ = EcPQ/(physicoChemical::k.value()*cloud_.constProps(typeIdQ).thetaV());

        idQ = cloud_.constProps(typeIdQ).thetaD()/cloud_.constProps(typeIdQ).thetaV();

        deltaDissoIQ = imaxQ-idQ;
        
        relax_ = true;

//         scalar P_exch = 0.0;

        if(deltaDissoIQ > 0)
        {
            // DISSOCIATION //

            nTotDissociationReactions_++;
            
            if (allowSplitting_)
            {
                relax_ = false; 
                
                scalar heatOfReactionDissJoules = heatOfReactionDiss_*physicoChemical::k.value();
                
                translationalEnergy = translationalEnergy + heatOfReactionDissJoules + EVibQ;
                
                //Particle Q is dissociating, so no vibrational redistribution to it
                
                //Particle P is an atom, so no redistribution required for it
                
//                 if (0.2 > cloud_.rndGen().scalar01())
//                 {
//                     scalar EcQ = translationalEnergy + ERotQ;
//                     
//                     scalar energyRatio = 0.0;
//                     
//                     if(rotationalDofQ == 2.0)
//                     {
//                         energyRatio = 1.0 - pow(cloud_.rndGen().scalar01(),(1.0/ChiB));
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

                label q1 = 0;
                label q2 = 1;

                const label& typeId1 = dissociationProductIds_[q1];
                const label& typeId2 = dissociationProductIds_[q2];
                
                //Mass of Product one and two
                scalar mP1 = cloud_.constProps(typeId1).mass();
                scalar mP2 = cloud_.constProps(typeId2).mass();
                
                scalar mRatoms = mP1*mP2/(mP1 + mP2);

//                 scalar cRatoms = sqrt(2.0*Ecplx/mRatoms);
                
                scalar cRatoms = sqrt(2.0*ERotQ/mRatoms);

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
                q.ERot() = 0.0;
                q.vibLevel() = 0;
                
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
                    0,
                    classificationQ
                );
            }
        }
        
        if(deltaDissoIQ <= 0)
        {
            scalar P_exch = 0.0;
            
            TColl = (translationalEnergy/(physicoChemical::k.value()))/(2.5 - omegaPQ);// Eq.(10) of Bird QK paper.
                
            scalar heatOfReactionExchJoules = heatOfReactionExch_*physicoChemical::k.value();

            scalar activationEnergy = activationEnergy_ + (aCoeff_*pow((TColl/273.0) , bCoeff_) * fabs(heatOfReactionExchJoules));
        
            if(EcPQ > activationEnergy) // i.e. exchange reaction can possibly occur.
            {
                scalar summation = 0.0; // declare "summation" term.

                if(activationEnergy < cloud_.constProps(typeIdQ).thetaV()*physicoChemical::k.value())
                {
                    summation = 1.0; // this refers to the first sentence in Bird's QK paper after Eq.(12).
                }
                else
                {
                    iaQ = EcPQ / (physicoChemical::k.value()*cloud_.constProps(typeIdQ).thetaV());

                    for(label i = 0 ; i <= iaQ ; i++)
                    {
                        summation += pow((1.0 - ((i*physicoChemical::k.value()*cloud_.constProps(typeIdQ).thetaV())/EcPQ)),(1.5 - omegaPQ));
                    }
                }
            
                P_exch = pow((1.0 - (activationEnergy/EcPQ)),(1.5 - omegaPQ))/summation;// now based on modified activation energy.
            }
        
            if(P_exch > cloud_.rndGen().scalar01())
            {   
                // EXCHANGE REACTION //

                nTotExchangeReactions_++;
                
                if(allowSplitting_)
                {   
                    relax_ = false;
                    
                    //center of mass velocity (pre-exchange)
                    vector Ucm = (mP*UP + mQ*UQ)/(mP + mQ);

                    const label& typeIdMol = exchangeProductIds_[0];
                    const label& typeIdAtom = exchangeProductIds_[1];

                    // change species properties

                    scalar mQExch = cloud_.constProps(typeIdAtom).mass();
                    scalar mPExch = cloud_.constProps(typeIdMol).mass();

                    scalar mRExch = mPExch*mQExch/(mPExch + mQExch);
                    
                    translationalEnergy = translationalEnergy + ERotQ + EVibQ + heatOfReactionExchJoules;
                    
                    p.typeId() = typeIdMol; // p is originally the atom, becomes the molecule
                    q.typeId() = typeIdAtom; // q is orinally the molecule, becomes the atom
                                        
                    scalar EVibP = 0.0;
                    
                    scalar EcP = translationalEnergy + EVibP;
                    label iMaxP = (EcP /(physicoChemical::k.value()*thetaV));

                    if(iMaxP > SMALL)
                    {
                        // - Bird equations 5.42 and 11.34 gave this denominator
                        scalar TCollP = (iMaxP*thetaV) / (3.5 - omegaPQExch); 
                        
                        scalar pow1 = pow((thetaD/TCollP),0.333) - 1.0;

                        scalar pow2 = pow ((thetaD/refTempZv),0.333) -1.0;
                        
                        // - vibrational collision number (equation 2, Bird 2010)
                        scalar ZvP1 = pow((thetaD/TCollP),omegaPQExch); 
                        
                        scalar ZvP2 = pow(Zref*(pow((thetaD/refTempZv),(-1.0*omegaPQExch))),(pow1/pow2));
                        
                        scalar ZvP = ZvP1*ZvP2;
//                         scalar ZvP = 50.0;
                            
                        scalar inverseVibrationalCollisionNumberP = 1.0/ZvP;
                    
                        if(inverseVibrationalCollisionNumberP > cloud_.rndGen().scalar01())
                        {
                            label iDashP = 0; // post-collision quantum number
                            scalar func = 0.0;

                            do // acceptance - rejection 
                            {
                                iDashP = cloud_.rndGen().integer(0,iMaxP);
                                EVibP = iDashP*physicoChemical::k.value()*thetaV;
                                func = pow((1.0 - (EVibP / EcP)),(1.5 - omegaPQExch));
                        
                            } while( !(func > cloud_.rndGen().scalar01()) );
                    
                            translationalEnergy = EcP - EVibP;
                        }
                    }
                    
                    scalar ERotP = 0.0;
                    
                    if (0.2 > cloud_.rndGen().scalar01())
                    {
                        scalar EcP = translationalEnergy + ERotP;
                        
                        scalar energyRatio = 0.0;
                        
                        if(rotationalDofExch == 2.0)
                        {
                            energyRatio = 1.0 - pow(cloud_.rndGen().scalar01(),(1.0/ChiBExch));
                        }
                        else
                        {
                            scalar ChiA = 0.5*rotationalDofExch;
                            
                            energyRatio = cloud_.energyRatio(ChiA, ChiBExch);
                        }

                        ERotP = energyRatio*EcP;
                    
                        translationalEnergy = EcP - ERotP;
                    }
                    
                    scalar relVelExchMol = sqrt((2.0*translationalEnergy)/mRExch);

                    // Variable Hard Sphere collision part for collision of molecules
            
                    scalar cosTheta = 2.0*cloud_.rndGen().scalar01() - 1.0;
                
                    scalar sinTheta = sqrt(1.0 - cosTheta*cosTheta);
                
                    scalar phi = twoPi*cloud_.rndGen().scalar01();
                
                    vector postCollisionRelU =
                        relVelExchMol
                        *vector
                            (
                                cosTheta,
                                sinTheta*cos(phi),
                                sinTheta*sin(phi)
                            );
                            
                    UP = Ucm + (postCollisionRelU*mQExch/(mPExch + mQExch)); // P changes from atom to mol.
                    UQ = Ucm - (postCollisionRelU*mPExch/(mPExch + mQExch)); // Q changes from mol to atom. 
                    
                    p.vibLevel() = EVibP/(cloud_.constProps(typeIdP).thetaV()*physicoChemical::k.value());
                    p.ERot() = ERotP;
                    p.U() = UP;
                    
                    q.U() = UQ;
                    q.ERot() = 0.0; // remove q's internal energy as it's now an atom
                    q.vibLevel() = 0.0;
                }
            }
        }
    }
}

void mixedTypeIIDissociationReverseExchange::reactExchangeMolecule
(
    dsmcParcel& p,
    label newTypeId,
    const label& newEVibLevel,
    const scalar& newERot,
    const vector& newU
)
{
    p.vibLevel() = newEVibLevel;
    p.ERot() = newERot;
    p.U() = newU;
    p.typeId() = newTypeId;
}

void mixedTypeIIDissociationReverseExchange::reactExchangeAtom
(
    dsmcParcel& p,
    label newTypeId,
    const vector& newU
)
{
    p.U() = newU;
    p.typeId() = newTypeId;
}

void  mixedTypeIIDissociationReverseExchange::outputResults(const label& counterIndex)
{
    if(writeRatesToTerminal_ == true)
    {
        
//         if(counterIndex == 1)
//         {
            const List< DynamicList<dsmcParcel*> >& cellOccupancy
                = cloud_.cellOccupancy();

            forAll(cellOccupancy, c)
            {
                volume_ += mesh_.cellVolumes()[c];
            }

            //- Parallel communication
//             if(Pstream::parRun())
//             {
//                 reduce(cellVolume_, sumOp<scalar>());
//             }
//         }
        
        // measure density 
//         if(counterIndex == 1)
//         {
//             const List< DynamicList<dsmcParcel*> >& cellOccupancy
//                 = cloud_.cellOccupancy();

            List<label> mols;
	    mols.append(0); mols.append(0);
            scalar volume = volume_;
            label nTotExchangeReactions = nTotExchangeReactions_;
            label nTotDissociationReactions = nTotDissociationReactions_;

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
                reduce(nTotExchangeReactions, sumOp<label>());
                reduce(nTotDissociationReactions, sumOp<label>());
            }

            numberDensities_[0] = (mols[0]*cloud().nParticle())/volume;
            numberDensities_[1] = (mols[1]*cloud().nParticle())/volume;
//         }

        const scalar& deltaT = mesh_.time().deltaT().value();

        word reactantMolA = cloud_.typeIdList()[reactantIds_[0]];
        word reactantMolB = cloud_.typeIdList()[reactantIds_[1]];

        word exchangeProductMolA = cloud_.typeIdList()[exchangeProductIds_[0]];
        word exchangeProductMolB = cloud_.typeIdList()[exchangeProductIds_[1]];

        word dissociationProductMolA = cloud_.typeIdList()[dissociationProductIds_[0]];
        word dissociationProductMolB = cloud_.typeIdList()[dissociationProductIds_[1]];

        if((numberDensities_[0] > 0.0) && (numberDensities_[1] > 0.0))
        {   

            scalar reactionRate1 = 0.0;
            scalar reactionRate2 = 0.0;

            reactionRate1 =
            (
                nTotExchangeReactions
                * cloud_.nParticle()
            )/(counterIndex*deltaT*numberDensities_[0]* numberDensities_[1]*volume);

            reactionRate2 =
            (
                nTotDissociationReactions
                * cloud_.nParticle()
            )/(counterIndex*deltaT*numberDensities_[0]* numberDensities_[1]*volume);

            Info << "Reverse exchange reaction " 
                <<  reactantMolA << " + " << reactantMolB 
                << " --> " << exchangeProductMolA << " + " << exchangeProductMolB
                << ", reaction rate = " << reactionRate1 << nl
                << "Type II dissociation reaction " 
                <<  reactantMolA << " + " << reactantMolB << " --> " 
                << dissociationProductMolA << " + " << dissociationProductMolB 
                <<  " + " << reactantMolB 
                << ", reaction rate = " << reactionRate2 
                << endl;

        }
    }
    else
    {
        word reactantMolA = cloud_.typeIdList()[reactantIds_[0]];
        word reactantMolB = cloud_.typeIdList()[reactantIds_[1]];

        word exchangeProductMolA = cloud_.typeIdList()[exchangeProductIds_[0]];
        word exchangeProductMolB = cloud_.typeIdList()[exchangeProductIds_[1]];

        word dissociationProductMolA = cloud_.typeIdList()[dissociationProductIds_[0]];
        word dissociationProductMolB = cloud_.typeIdList()[dissociationProductIds_[1]];
        
        label nTotExchangeReactions = nTotExchangeReactions_;
        label nTotDissociationReactions = nTotDissociationReactions_;
        
        //- Parallel communication
        if(Pstream::parRun())
        {
            reduce(nTotExchangeReactions, sumOp<label>());
            reduce(nTotDissociationReactions, sumOp<label>());
        }
        
        if(nTotExchangeReactions > VSMALL)
        {
            Info << "Reverse exchange reaction " 
                    <<  reactantMolA << " + " << reactantMolB << " --> " 
                    << exchangeProductMolA << " + " << exchangeProductMolB 
                    << " is active." << endl;
        }
        
        if(nTotDissociationReactions > VSMALL)
        {
            Info  << "Type II dissociation reaction " 
                <<  reactantMolA << " + " << reactantMolB << " --> " 
                << dissociationProductMolA << " + " << dissociationProductMolB 
                <<  " + " << reactantMolB 
                << " is active." << endl;
        }
    }
    
    nReactionsPerTimeStep_ = 0.0;
}


const bool& mixedTypeIIDissociationReverseExchange::relax() const
{
    return relax_;
}

}
// End namespace Foam

// ************************************************************************* //
