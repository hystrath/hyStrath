/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2007 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.
f
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

#include "mixedTypeIIDissociationForwardExchange.H"
#include "addToRunTimeSelectionTable.H"
#include "fvc.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(mixedTypeIIDissociationForwardExchange, 0);

addToRunTimeSelectionTable(dsmcReaction, mixedTypeIIDissociationForwardExchange, dictionary);



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
mixedTypeIIDissociationForwardExchange::mixedTypeIIDissociationForwardExchange
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
    cellVolume_(0.0),
    numberDensities_(2, 0.0)
{

}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

mixedTypeIIDissociationForwardExchange::~mixedTypeIIDissociationForwardExchange()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void mixedTypeIIDissociationForwardExchange::initialConfiguration()
{
    setProperties();
}

void mixedTypeIIDissociationForwardExchange::setProperties()
{
    // reading in reactants


    const List<word> reactantMolecules (propsDict_.lookup("reactantMolecules"));

    if(reactantMolecules.size() != 2)
    {
        FatalErrorIn("mixedTypeIIDissociationForwardExchange::setProperties()")
            << "There should be two or more reactants, instead of " 
            << reactantMolecules.size() << nl 
            << exit(FatalError);
    }
    
    if(reactantMolecules[0] == reactantMolecules[1])
    {
        FatalErrorIn("mixedTypeIIDissociationForwardExchange::setProperties()")
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
            FatalErrorIn("mixedTypeIIDissociationForwardExchange::setProperties()")
                << "Cannot find type id: " << reactantMolecules[r] << nl 
                << exit(FatalError);
        }
    }
    
    // check that the first reactants is a 'MOLECULE' 

    const scalar& rDof1 = cloud_.constProps(reactantIds_[0]).rotationalDegreesOfFreedom();

    if(rDof1 < 1)
    {
	FatalErrorIn("mixedTypeIIDissociationForwardExchange::setProperties()")
	    << "Reactant 1 must be a molecule (not an atom): " << reactantMolecules[0] 
	    << nl 
	    << exit(FatalError);
    }
    
    // check that the second reactant is an 'ATOM' 

    const scalar& rDof2 = cloud_.constProps(reactantIds_[1]).rotationalDegreesOfFreedom();

    if(rDof2 > 0)
    {
	FatalErrorIn("mixedTypeIIDissociationForwardExchange::setProperties()")
	    << "Reactant 2 must be an atom (not a molecule): " << reactantMolecules[1] 
	    << nl 
	    << exit(FatalError);
    }

    // reading in products

    const List<word> exchangeProductMolecules (propsDict_.lookup("productsOfExchangeReaction"));

    if(exchangeProductMolecules.size() != 2)
    {
        FatalErrorIn("mixedTypeIIDissociationForwardExchange::setProperties()")
            << "There should be two mixedTypeIIDissociationForwardExchange reaction products, instead of " 
            << reactantMolecules.size() << nl 
            << exit(FatalError);
    }
    
    if(exchangeProductMolecules[0] == exchangeProductMolecules[1])
    {
        FatalErrorIn("mixedTypeIIDissociationForwardExchange::setProperties()")
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
            FatalErrorIn("mixedTypeIIDissociationForwardExchange::setProperties()")
                << "Cannot find type id: " << exchangeProductMolecules[r] << nl 
                << exit(FatalError);
        }
    }
    
    const List<word> dissociationProductMolecules (propsDict_.lookup("productsOfDissociatedMolecule"));
    
    if(dissociationProductMolecules.size() != 2)
    {
        FatalErrorIn("mixedTypeIIDissociationForwardExchange::setProperties()")
            << "There should be two mixedTypeIIDissociationForwardExchange reaction products, instead of " 
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
            FatalErrorIn("mixedTypeIIDissociationForwardExchange::setProperties()")
                << "Cannot find type id: " << dissociationProductMolecules[r] << nl 
                << exit(FatalError);
        }
    }
    
    //check that first exchange product is a 'MOLECULE' (not an 'ATOM')
    
    const scalar& rDof3 = cloud_.constProps(exchangeProductIds_[0]).rotationalDegreesOfFreedom();

    if(rDof3 < 0)
    {
        FatalErrorIn("mixedTypeIIDissociationForwardExchange::setProperties()")
            << "First product of the exchange reaction must be a molecule (not an atom): " << exchangeProductMolecules[0] 
            << nl 
            << exit(FatalError);
    }
    
    //check that second exchange product is an 'ATOM' (not a 'MOLECULE')
    
    const scalar& rDof4 = cloud_.constProps(exchangeProductIds_[1]).rotationalDegreesOfFreedom();

    if(rDof4 > 0)
    {
        FatalErrorIn("mixedTypeIIDissociationForwardExchange::setProperties()")
            << "Second product of the exchange reaction must be an atom (not a molecule): " << exchangeProductMolecules[1] 
            << nl 
            << exit(FatalError);
    }
    
        const scalar& rDof5= cloud_.constProps(dissociationProductIds_[0]).rotationalDegreesOfFreedom();

    if(rDof5 > 0)
    {
        FatalErrorIn("mixedTypeIIDissociationForwardExchange::setProperties()")
            << "First product of the dissociation reaction must be an atom (not a molecule): " << dissociationProductMolecules[0] 
            << nl 
            << exit(FatalError);
    }
    
    //check that second exchange product is an 'ATOM' (not a 'MOLECULE')
    
    const scalar& rDof6 = cloud_.constProps(dissociationProductIds_[1]).rotationalDegreesOfFreedom();

    if(rDof6 > 0)
    {
        FatalErrorIn("mixedTypeIIDissociationForwardExchange::setProperties()")
            << "Second product of the exchange reaction must be an atom (not a molecule): " << dissociationProductMolecules[1] 
            << nl 
            << exit(FatalError);
    }

    activationEnergy_ *= physicoChemical::k.value();
}

bool mixedTypeIIDissociationForwardExchange::tryReactMolecules(const label& typeIdP, const label& typeIdQ)
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

void mixedTypeIIDissociationForwardExchange::reaction
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


void mixedTypeIIDissociationForwardExchange::reaction
(
    dsmcParcel& p,
    dsmcParcel& q
)
{

    label typeIdP = p.typeId();
    label typeIdQ = q.typeId();

    if(typeIdP == reactantIds_[0] && typeIdQ == reactantIds_[1]) // AB + C ---> A + B + C (diss type II) or AB + C ---> AC + B (forward exchange)
    {
        vector UP = p.U();
        vector UQ = q.U();
        scalar ERotP = p.ERot();
        scalar EVibP = p.EVib();

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

        scalar P_diss = 0.0;
        scalar P_exch = 0.0;

        if(deltaDissoIP > 0)
        {
            P_diss = 1.0;
        }

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


        if( (2.0*cloud_.rndGen().scalar01()) < (P_diss + P_exch))
//         if( (cloud_.rndGen().scalar01()) < (P_diss + P_exch))
        {
            scalar P_react = (P_exch / (P_exch + P_diss));

            if(P_react > (0.5*cloud_.rndGen().scalar01()))            
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

                    scalar mP = cloud_.constProps(typeIdAtom).mass();

                    scalar mQ = cloud_.constProps(typeIdMol).mass();

                    scalar thetaVQ = cloud_.constProps(typeIdMol).thetaV();

                    scalar mR = mP*mQ/(mP + mQ);

                    scalar omegaPQ =
                    0.5
                    *(
                        cloud_.constProps(typeIdAtom).omega()
                        + cloud_.constProps(typeIdMol).omega()
                    );

                    // Exchange internal energy of the parcels

                    scalar Erot2 = ERotP;
                    scalar Evib2 = EVibP;

                    //////////////////////////////////////////////////////////////////////////////////
                    // Determine the pre-collision total energy. Note that for this endothermic     //
                    // exchange reaction the heatOfReactionJoules is subtracted from the total. //
                    // heatOfReactionJoules is negative in the chemReactDict.                       //
                    //////////////////////////////////////////////////////////////////////////////////

                    scalar heatOfReactionJoules = heatOfReactionExch_*physicoChemical::k.value();

                    scalar EcTot = fabs(translationalEnergy + Erot2 + Evib2 + heatOfReactionJoules);

                    scalar func = 0.0;
                    scalar rel = 0.0;
                    scalar psiV = 0.0;

                    ///////////////////////////////////////
                    // acceptance-rejection method (ARM) //
                    ///////////////////////////////////////

                    scalar DOFrot1 = cloud_.constProps(exchangeProductIds_[0]).rotationalDegreesOfFreedom();
                    scalar DOFrot2 = cloud_.constProps(exchangeProductIds_[1]).rotationalDegreesOfFreedom();

                    scalar DOFm = 5.0 - (2.0*omegaPQ);
                    scalar DOFsum = DOFm + DOFrot1 + DOFrot2;

                    scalar expo = 0.5*DOFsum -1.0;

                    rel = EcTot / (physicoChemical::k.value()*thetaVQ);
                    
                    do
                    {
                        psiV = label(((label(rel)) +1)*cloud_.rndGen().scalar01())/rel;
                        func = pow((1.0-psiV),expo);
                    } while (func < cloud_.rndGen().scalar01());

                    Evib2 = psiV* EcTot;

                    EcTot -= Evib2;

                    // Distribute EcTot over translational and all rotational modes of P and Q based on MONACO code //

                    DOFm = 1.0*(5.0 - (2.0 * omegaPQ)); // part of MONACO code

                    scalar DOFtot = DOFm + DOFrot1 + DOFrot2;

                    scalar rPSIm = 0.0;

                    scalar h1 = 0.5*DOFtot - 2.0;
                    scalar h2 = 0.5*DOFm - 1.0 + SMALL;
                    scalar h3 = 0.5*(DOFtot-DOFm)-1.0 + SMALL;

                    scalar prob = 0.0;

                    do
                    {
                        rPSIm = cloud_.rndGen().scalar01();
                        prob = pow(h1,h1)/(pow(h2,h2)*pow(h3,h3))*pow(rPSIm,h2)*pow(1.0-rPSIm,h3);

                    } while (prob < cloud_.rndGen().scalar01());

                    if (DOFm == DOFtot) rPSIm = 1.0;

                    if (DOFm == 2.0 && DOFtot == 4.0) rPSIm = cloud_.rndGen().scalar01();

                    if (DOFtot < 4.0 ) rPSIm = (DOFm/DOFtot);

                    scalar Erel = EcTot*rPSIm;
                    scalar Erot = EcTot - Erel;

                    Erot2 = Erot;

                    scalar relVelExchMol = pow(2*Erel/mR,0.5);

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
            
                    UP = Ucm + (postCollisionRelU*mQ/(mP + mQ)); // P changes from mol to atom.
                    UQ = Ucm - (postCollisionRelU*mP/(mP + mQ)); // Q changes from atom to mol. 

                    if (cloud_.constProps(q.typeId()).rotationalDegreesOfFreedom() > VSMALL)
                    { 
                        q.EVib() = Evib2;
                        q.ERot() = Erot2;
                        q.U() = UQ;
                        q.typeId() = typeIdMol;
                        p.U() = UP;
                        p.typeId() = typeIdAtom;
                    }
                    else
                    { 
                        p.EVib() = Evib2;
                        p.ERot() = Erot2;
                        p.U() = UQ;
                        p.typeId() = typeIdMol;
                        q.U() = UP;
                        q.typeId() = typeIdAtom;
                    }
                }
            }
            else
            {
                // DISSOCIATION //

                nTotDissociationReactions_++;
                
                if (allowSplitting_)
                {
                    relax_ = false;
                    
                    //////////////////////////////////////////////////////////////////////////////////
                    // Determine the pre-collision total energy. Note that for this endothermic     //
                    // dissociation reaction the heatOfReactionJoules is subtracted from the total. //
                    // heatOfReactionJoules is negative in the chemReactDict.                       //
                    //////////////////////////////////////////////////////////////////////////////////

                    scalar heatOfReactionJoules = heatOfReactionDiss_*physicoChemical::k.value();
                    scalar EcTot = fabs(translationalEnergy + ERotP + EVibP + heatOfReactionJoules);

                    // Distribute EcTot over translational and all rotational mode of P based on MONACO code //


                    scalar DOFm = 2.0*(5.0 - (2.0 * omegaPQ)); // 2.0* at start is in 1st part of MONACO code

                    scalar DOFrot1 = cloud_.constProps(reactantIds_[0]).rotationalDegreesOfFreedom();
                    scalar DOFrot2 = cloud_.constProps(reactantIds_[1]).rotationalDegreesOfFreedom();

                    scalar DOFtot = DOFm + DOFrot1 + DOFrot2;

                    scalar h1 = 0.5*DOFtot - 2.0;
                    scalar h2 = 0.5*DOFm - 1.0 + 1.0e-5;
                    scalar h3 = 0.5*(DOFtot-DOFm)-1.0 + 1.0e-5;

                    scalar rPSIm = 0.0;
                    scalar prob = 0.0;

                    do
                    {
                        rPSIm = cloud_.rndGen().scalar01();
                        prob = pow(h1,h1)/(pow(h2,h2)*pow(h3,h3))*pow(rPSIm,h2)*pow(1.0-rPSIm,h3);

                    } while (prob < cloud_.rndGen().scalar01());

                    if (DOFm == DOFtot) rPSIm = 1.0;

                    if (DOFm == 2.0 && DOFtot == 4.0) rPSIm = cloud_.rndGen().scalar01();

                    if (DOFtot < 4.0 ) rPSIm = (DOFm/DOFtot);

                    scalar Erel = EcTot*rPSIm;
                    scalar Erot = EcTot - Erel;

                    scalar Ecplx = Erot;

                    scalar relVelNonDissoMol = pow(2*(Erel)/mR,0.5);

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

                    scalar cRatoms = sqrt(2.0* Ecplx/mRatoms);

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
                    label cellI = p.cell();
                    vector position = p.position();
                    label tetFaceI = p.tetFace();
                    label tetPtI = p.tetPt();
                    
                    p.typeId() = typeId1;
                    p.U() = uP1;
                    p.EVib() = 0.0;
                    p.ERot() = 0.0;
                    
                    // insert new product 2
                    cloud_.addNewParcel
                    (
                        position,
                        uP2,
                        0.0,
                        0.0,
                        cellI,
                        tetFaceI,
                        tetPtI,
                        typeId2,
                        0
                    );
                }
            }
        }
        else
        {
            relax_ = true; // no reaction - use LB
        }

    } // end of typeIdP
     
    if(typeIdP == reactantIds_[1] && typeIdQ == reactantIds_[0])// C + AB ---> A + B + C (diss type II) or C + AB ---> AC + B (forward exchange)
    {
        vector UP = p.U();
        vector UQ = q.U();
        scalar ERotQ = q.ERot();
        scalar EVibQ = q.EVib();

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

        scalar P_diss = 0.0;
        scalar P_exch = 0.0;

        if(deltaDissoIQ > 0)
        {
            P_diss = 1.0;
        }

        // Next, calculate exchange probability.

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


        if ( (2.0*cloud_.rndGen().scalar01()) < (P_diss + P_exch))
//         if ( (cloud_.rndGen().scalar01()) < (P_diss + P_exch))
        {
            scalar P_react = (P_exch / (P_exch + P_diss));

            if (P_react > (0.5*cloud_.rndGen().scalar01()))
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

                    scalar mP = cloud_.constProps(typeIdMol).mass();

                    scalar mQ = cloud_.constProps(typeIdAtom).mass();

                    scalar thetaVP = cloud_.constProps(typeIdMol).thetaV();

                    scalar mR = mP*mQ/(mP + mQ);

                    scalar omegaPQ =
                    0.5
                    *(
                        cloud_.constProps(typeIdAtom).omega()
                        + cloud_.constProps(typeIdMol).omega()
                    );

                    // Exchange internal energy of the parcels

                    scalar Erot2 = ERotQ;
                    scalar Evib2 = EVibQ;

                    //////////////////////////////////////////////////////////////////////////////////
                    // Determine the pre-collision total energy. Note that for this endothermic     //
                    // exchange reaction the heatOfReactionJoules is subtracted from the total. //
                    // heatOfReactionJoules is negative in the chemReactDict.                       //
                    //////////////////////////////////////////////////////////////////////////////////

                    scalar heatOfReactionJoules = heatOfReactionExch_ * physicoChemical::k.value();

                    scalar EcTot = fabs(translationalEnergy + Erot2 + Evib2 + heatOfReactionJoules);

                    scalar func = 0.0;
                    scalar rel = 0.0;
                    scalar psiV = 0.0;

                    ///////////////////////////////////////
                    // acceptance-rejection method (ARM) //
                    ///////////////////////////////////////

                    scalar DOFrot1 = cloud_.constProps(exchangeProductIds_[0]).rotationalDegreesOfFreedom();
                    scalar DOFrot2 = cloud_.constProps(exchangeProductIds_[1]).rotationalDegreesOfFreedom();

                    scalar DOFm = 5.0 - (2.0*omegaPQ);
                    scalar DOFsum = DOFm + DOFrot1 + DOFrot2;

                    scalar expo = 0.5*DOFsum -1.0;

                    rel = EcTot / (physicoChemical::k.value()*thetaVP);
                    
                    do
                    {
                        psiV = label(((label(rel)) +1)*cloud_.rndGen().scalar01())/rel;
                        func = pow((1.0-psiV),expo);
                    } while (func < cloud_.rndGen().scalar01());

                    Evib2 = psiV* EcTot;

                    EcTot -= Evib2;


                    // Distribute EcTot over translational and all rotational modes of P and Q based on MONACO code //

                    DOFm = 1.0*(5.0 - (2.0 * omegaPQ)); // part of MONACO code

                    scalar DOFtot = DOFm + DOFrot1 + DOFrot2;

                    scalar rPSIm = 0.0;

                    scalar h1 = 0.5*DOFtot - 2.0;
                    scalar h2 = 0.5*DOFm - 1.0 + SMALL;
                    scalar h3 = 0.5*(DOFtot-DOFm)-1.0 + SMALL;

                    scalar prob = 0.0;

                    do
                    {
                        rPSIm = cloud_.rndGen().scalar01();
                        prob = pow(h1,h1)/(pow(h2,h2)*pow(h3,h3))*pow(rPSIm,h2)*pow(1.0-rPSIm,h3);

                    } while (prob < cloud_.rndGen().scalar01());

                    if (DOFm == DOFtot) rPSIm = 1.0;

                    if (DOFm == 2.0 && DOFtot == 4.0) rPSIm = cloud_.rndGen().scalar01();

                    if (DOFtot < 4.0 ) rPSIm = (DOFm/DOFtot);

                    scalar Erel = EcTot*rPSIm;
                    scalar Erot = EcTot - Erel;

                    Erot2 = Erot;

                    scalar relVelExchMol = pow(2*Erel/mR,0.5);

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
            
                    UP = Ucm + (postCollisionRelU*mQ/(mP + mQ)); // P changes from atom to mol.
                    UQ = Ucm - (postCollisionRelU*mP/(mP + mQ)); // Q changes from mol to atom. 

                    if (cloud_.constProps(p.typeId()).rotationalDegreesOfFreedom() > VSMALL)
                    {
                        p.EVib() = Evib2;
                        p.ERot() = Erot2;
                        p.U() = UP;
                        p.typeId() = typeIdMol;
                        q.U() = UQ;
                        q.typeId() = typeIdAtom;
                    }
                    else   
                    {
                        q.EVib() = Evib2;
                        q.ERot() = Erot2;
                        q.U() = UP;
                        q.typeId() = typeIdMol;
                        p.U() = UQ;
                        p.typeId() = typeIdAtom;
                    }
                }
            }
            else
            {
                // DISSOCIATION //

                nTotDissociationReactions_++;
                
                if (allowSplitting_)
                {
                    relax_ = false;

                    //////////////////////////////////////////////////////////////////////////////////
                    // Determine the pre-collision total energy. Note that for this endothermic     //
                    // dissociation reaction the heatOfReactionJoules is subtracted from the total. //
                    // heatOfReactionJoules is negative in the chemReactDict.                       //
                    //////////////////////////////////////////////////////////////////////////////////

                    scalar heatOfReactionJoules = heatOfReactionDiss_*physicoChemical::k.value();

                    scalar EcTot = fabs(translationalEnergy + ERotQ + EVibQ + heatOfReactionJoules);

                    // Distribute EcTot over translational and all rotational mode of P based on MONACO code //

                    scalar DOFm = 2.0*(5.0 - (2.0 * omegaPQ)); // 2.0* at start is in 1st part of MONACO code now reduced

                    scalar DOFrot1 = cloud_.constProps(reactantIds_[0]).rotationalDegreesOfFreedom();
                    scalar DOFrot2 = cloud_.constProps(reactantIds_[1]).rotationalDegreesOfFreedom();

                    scalar DOFtot = DOFm + DOFrot1 + DOFrot2;

                    scalar h1 = 0.5*DOFtot - 2.0;
                    scalar h2 = 0.5*DOFm - 1.0 + SMALL;
                    scalar h3 = 0.5*(DOFtot-DOFm)-1.0 + SMALL;

                    scalar rPSIm = 0.0;
                    scalar prob = 0.0;

                    do
                    {
                        rPSIm = cloud_.rndGen().scalar01();
                        prob = pow(h1,h1)/(pow(h2,h2)*pow(h3,h3))*pow(rPSIm,h2)*pow(1.0-rPSIm,h3);

                    } while (prob < cloud_.rndGen().scalar01());

                    if (DOFm == DOFtot) rPSIm = 1.0;

                    if (DOFm == 2.0 && DOFtot == 4.0) rPSIm = cloud_.rndGen().scalar01();

                    if (DOFtot < 4.0 ) rPSIm = (DOFm/DOFtot);

                    scalar Erel = EcTot*rPSIm;
                    scalar Erot = EcTot - Erel;

                    scalar Ecplx = Erot;

                    scalar relVelNonDissoMol = pow(2*Erel/mR,0.5);

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


                    scalar cRatoms = sqrt(2.0* Ecplx/mR);

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

                    label q1 = 0;
                    label q2 = 1;

                    const label& typeId1 = dissociationProductIds_[q1];
                    const label& typeId2 = dissociationProductIds_[q2];

                    //Mass of Product one and two
                    scalar mP1 = cloud_.constProps(typeId1).mass();
                    scalar mP2 = cloud_.constProps(typeId2).mass();

                    vector uP1 = UQ + postCollisionRelU2*mP2/(mP1 + mP2);
                    vector uP2 = UQ - postCollisionRelU2*mP1/(mP1 + mP2);

                    // P remains NON-DISSOCIATED (it is an atom).
                    p.U() = UP;

                    // Molecule Q will dissociate into 2 atoms.
                    label cellI = p.cell();
                    vector position = p.position();
                    label tetFaceI = p.tetFace();
                    label tetPtI = p.tetPt();
                    
//                     cloud_.deleteParticle(q);
//             
//                     // insert new product 1
//                     cloud_.addNewParcel
//                     (
//                         position,
//                         uP1,
//                         0.0,
//                         0.0,
//                         cellI,
//                         tetFaceI,
//                         tetPtI,
//                         typeId1,
//                         0
//                     );
                    
                    q.typeId() = typeId1;
                    q.U() = uP1;
                    q.EVib() = 0.0;
                    q.ERot() = 0.0;
                    
                    // insert new product 2
                    cloud_.addNewParcel
                    (
                        position,
                        uP2,
                        0.0,
                        0.0,
                        cellI,
                        tetFaceI,
                        tetPtI,
                        typeId2,
                        0
                    );
                }
            } // end of allowSplitting_
        }
        else
        {
            relax_ = true; // no reaction - use LB
        }

	} // end of typeIdP

} // end of reaction class

void mixedTypeIIDissociationForwardExchange::reactExchangeMolecule
(
    dsmcParcel& p,
    label newTypeId,
    const scalar& newEVib,
    const scalar& newERot,
    const vector& newU
)
{
    p.EVib() = newEVib;
    p.ERot() = newERot;
    p.U() = newU;
    p.typeId() = newTypeId;
}

void mixedTypeIIDissociationForwardExchange::reactExchangeAtom
(
    dsmcParcel& p,
    label newTypeId,
    const vector& newU
)
{
    p.U() = newU;
    p.typeId() = newTypeId;
}

void  mixedTypeIIDissociationForwardExchange::outputResults(const label& counterIndex)
{
    if(writeRatesToTerminal_ == true)
    {           
        // measure density 
        if(counterIndex == 1)
        {
            const List< DynamicList<dsmcParcel*> >& cellOccupancy
                = cloud_.cellOccupancy();

            List<label> mols;
	mols.append(0);
	mols.append(0);

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

                cellVolume_ += mesh_.cellVolumes()[c];
            }

            //- Parallel communication
            if(Pstream::parRun())
            {
                reduce(cellVolume_, sumOp<scalar>());
                reduce(mols[0], sumOp<label>());
                reduce(mols[1], sumOp<label>());
            }

            numberDensities_[0] = (mols[0]*cloud().nParticle())/cellVolume_;
            numberDensities_[1] = (mols[1]*cloud().nParticle())/cellVolume_;
        }


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
                nTotExchangeReactions_
                * cloud_.nParticle()
            )/(counterIndex*deltaT*numberDensities_[0]* numberDensities_[1]*cellVolume_);

            reactionRate2 =
            (
                nTotDissociationReactions_
                * cloud_.nParticle()
                )/(counterIndex*deltaT*numberDensities_[0]* numberDensities_[1]*cellVolume_);

            Info << "Forward exchange reaction " 
                    <<  reactantMolA << " + " << reactantMolB << " <--> " 
                    << exchangeProductMolA << " + " << exchangeProductMolB 
                    << ", reaction rate = " << reactionRate1 
                    << nl
                    << "Type II dissociation reaction " 
                    <<  reactantMolA << " + " << reactantMolB << " --> " 
                    << dissociationProductMolA << " + " << dissociationProductMolB 
                    <<  " + " << reactantMolB 
                    << ", reaction rate = " << reactionRate2
                    << endl;
        }
    }

    nReactionsPerTimeStep_ = 0.0;

}


const bool& mixedTypeIIDissociationForwardExchange::relax() const
{
    return relax_;
}

}
// End namespace Foam

// ************************************************************************* //
