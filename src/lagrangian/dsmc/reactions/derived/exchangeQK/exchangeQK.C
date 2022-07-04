/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016-2021 hyStrath
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of hyStrath, a derivative work of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Description

\*---------------------------------------------------------------------------*/

#include "exchangeQK.H"
#include "addToRunTimeSelectionTable.H"
#include "fvc.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(exchangeQK, 0);

addToRunTimeSelectionTable(dsmcReaction, exchangeQK, dictionary);


// * * * * * * * * * * *  Protected Member functions * * * * * * * * * * * * //

void exchangeQK::setProperties()
{
    dsmcReaction::setProperties();

    if (reactantIds_.size() != 2)
    {
        //- There must be exactly 2 reactants
        FatalErrorIn("exchangeQK::setProperties()")
            << "For reaction named " << reactionName_ << nl
            << "There should be two reactants, instead of "
            << reactantIds_.size() << nl
            << exit(FatalError);
    }

    bool moleculeFound = false;
    bool atomFound = false;

    forAll(reactantIds_, r)
    {
        //- Check if this reactant is a molecule
        if (reactantTypes_[r] >= 20)
        {
            moleculeFound = true;
            posMolReactant_ = r;
        }
        //- Check if this reactant is an atom
        else if (reactantTypes_[r] == 10 or reactantTypes_[r] == 11)
        {
            atomFound = true;
        }
        else
        {
            FatalErrorIn("exchangeQK::setProperties()")
                << "For reaction named " << reactionName_ << nl
                << "Reactant " << cloud_.typeIdList()[reactantIds_[r]]
                << " is neither a molecule nor an atom" << nl
                << exit(FatalError);
        }
    }

    if (!moleculeFound)
    {
        FatalErrorIn("exchangeQK::setProperties()")
            << "For reaction named " << reactionName_ << nl
            << "None of the reactants is a molecule." << nl
            << exit(FatalError);
    }
    else if (!atomFound)
    {
        FatalErrorIn("exchangeQK::setProperties()")
            << "For reaction named " << reactionName_ << nl
            << "None of the reactants is an atom." << nl
            << exit(FatalError);
    }

    //- Reading in exchange products
    const wordList productsExchange(propsDict_.lookup("exchangeProducts"));

    if (productsExchange.size() != 2)
    {
        FatalErrorIn("exchangeQK::setProperties()")
            << "For reaction named " << reactionName_ << nl
            << "There should be two products, instead of "
            << productsExchange.size() << nl
            << exit(FatalError);
    }

    productIdsExchange_.setSize(productsExchange.size());

    moleculeFound = false;
    atomFound = false;

    forAll(productIdsExchange_, r)
    {
        const label productIndex =
            findIndex
            (
                cloud_.typeIdList(),
                productsExchange[r]
            );

        //- Check that products belong to the typeIdList as defined in
        //  constant/dsmcProperties
        if (productIndex == -1)
        {
            FatalErrorIn("exchangeQK::setProperties()")
                << "For reaction named " << reactionName_ << nl
                << "Cannot find type id: " << productsExchange[r] << nl
                << exit(FatalError);
        }

        //- Check if this product is a molecule
        if (cloud_.constProps(productIndex).type() >= 20)
        {
            moleculeFound = true;
            //- The molecule is set to be the first product
            productIdsExchange_[0] = productIndex;
        }
        //- Check if this product is an atom
        else if
        (
            cloud_.constProps(productIndex).type() == 10
         or cloud_.constProps(productIndex).type() == 11
        )
        {
            atomFound = true;
            //- The atom is set to be the second product
            productIdsExchange_[1] = productIndex;
        }
        else
        {
            FatalErrorIn("exchangeQK::setProperties()")
                << "For reaction named " << reactionName_ << nl
                << "Product " << cloud_.typeIdList()[productIndex]
                << " is neither a molecule nor an atom" << nl
                << exit(FatalError);
        }
    }

    if (!moleculeFound)
    {
        FatalErrorIn("exchangeQK::setProperties()")
            << "For reaction named " << reactionName_ << nl
            << "None of the products is a molecule." << nl
            << exit(FatalError);
    }
    else if (!atomFound)
    {
        FatalErrorIn("exchangeQK::setProperties()")
            << "For reaction named " << reactionName_ << nl
            << "None of the products is an atom." << nl
            << exit(FatalError);
    }
}


void exchangeQK::testExchange
(
    const dsmcParcel& p,
    const scalar translationalEnergy,
    const scalar omegaPQ,
    scalar& collisionEnergy,
    scalar& totalReactionProbability,
    scalar& reactionProbability
)
{
    const label typeIdP = p.typeId();
    
    const scalar chiB = 2.5 - omegaPQ;

    //- Collision temperature: Eq.(10) of Bird's QK paper.
    const scalar TColl = translationalEnergy/(physicoChemical::k.value()*chiB);

    const scalar aDash =
        aCoeff_
       *(
            pow(chiB, bCoeff_)*exp(lgamma(chiB))/exp(lgamma(chiB + bCoeff_))
        );

    scalar activationEnergy =
        (
            aDash*pow(TColl/273.0, bCoeff_)
           *fabs(heatOfReactionExchangeJoules_)
        );
        
    scalar summation = 1.0;    

    if (heatOfReactionExchangeJoules_ < 0.0)
    {
        //- forward (endothermic) exchange reaction
        activationEnergy -= heatOfReactionExchangeJoules_;
    }

    label m = 0;
    
    do
    {
        const label vibLevel_m = p.vibLevel()[m];
        const scalar kBByThetaVP = physicoChemical::k.value()*cloud_.constProps(typeIdP).thetaV_m(m);
        const scalar EVibP_m = cloud_.constProps(typeIdP).eVib_m(m, vibLevel_m);

        //- Total collision energy
        collisionEnergy = translationalEnergy + EVibP_m;

        //- Condition for the exchange reaction to possibly occur
        if (collisionEnergy > activationEnergy)
        {
            if (activationEnergy > kBByThetaVP)
            {
                summation = 0.0;
                
                const label iaP = collisionEnergy/kBByThetaVP;

                for(label i=0; i<=iaP; i++)
                {
                    summation +=
                        pow
                        (
                            1.0 - cloud_.constProps(typeIdP).eVib_m(m, i)/collisionEnergy,
                            1.5 - omegaPQ
                        );
                }
            }

            //- Based on modified activation energy
            reactionProbability =
                pow
                (
                    1.0 - activationEnergy/collisionEnergy,
                    1.5 - omegaPQ
                )
                /summation;

            // Condition to exit the do while loop
            m = p.vibLevel().size();
        }

        m += 1;

    } while (m < p.vibLevel().size());

    totalReactionProbability += reactionProbability;
}


void exchangeQK::exchange
(
    dsmcParcel& p,
    dsmcParcel& q,
    scalar collisionEnergy
)
{
    const label typeIdP = p.typeId();
    const label typeIdQ = q.typeId();

    nTotExchangeReactions_++;
    nExchangeReactionsPerTimeStep_++;

    if (allowSplitting_)
    {
        relax_ = false;

        vector UP = p.U();
        vector UQ = q.U();

        const scalar mP = cloud_.constProps(typeIdP).mass();
        const scalar mQ = cloud_.constProps(typeIdQ).mass();
        const scalar mR = mP*mQ/(mP + mQ);
        const scalar cRsqr = magSqr(UP - UQ);

        scalar translationalEnergy = 0.5*mR*cRsqr;

        //- Center of mass velocity (pre-exchange)
        const vector Ucm = (mP*UP + mQ*UQ)/(mP + mQ);

        const label typeIdMol = productIdsExchange_[0];
        const label typeIdAtom = productIdsExchange_[1];

        //- Change species properties
        const scalar mPExch = cloud_.constProps(typeIdAtom).mass();
        const scalar mQExch = cloud_.constProps(typeIdMol).mass();
        const scalar mRExch = mPExch*mQExch/(mPExch + mQExch);
        const scalar omegaExch =
            0.5
            *(
                  cloud_.constProps(typeIdAtom).omega()
                + cloud_.constProps(typeIdAtom).omega()
            );

        const scalar EVibP = cloud_.constProps(typeIdP).eVib_tot(p.vibLevel());
        const scalar EEleP = cloud_.constProps(typeIdP).electronicEnergyList()[p.ELevel()];
        const scalar EEleQ = cloud_.constProps(typeIdQ).electronicEnergyList()[q.ELevel()];

        translationalEnergy += p.ERot() + EVibP + EEleP + EEleQ
            + heatOfReactionExchangeJoules_;

        //- p is originally the molecule and becomes the atom
        p.typeId() = typeIdAtom;
        p.ERot() = 0.0;
        p.vibLevel().setSize
        (
            cloud_.constProps
            (
                typeIdAtom
            ).nVibrationalModes(),
            0
        );
        p.ELevel() = 0;

        //- q is originally the atom and becomes the molecule
        q.typeId() = typeIdMol;
        q.ERot() = 0.0;
        q.vibLevel().setSize
        (
            cloud_.constProps
            (
                typeIdMol
            ).nVibrationalModes(),
            0
        );
        q.ELevel() = 0;
        
        //- Energy redistribution for particle q
        cloud_.binaryCollision().redistribute
        (
            q, translationalEnergy, omegaExch, true
        );
        
        //- Post-collision velocities
        const scalar relVelExchMol = sqrt(2.0*translationalEnergy/mRExch);
        
        p.U() = Ucm;
        
        cloud_.binaryCollision().postReactionVelocities
        (
            typeIdAtom,
            typeIdMol,
            p.U(),
            q.U(),
            relVelExchMol
        );
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
exchangeQK::exchangeQK
(
    Time& t,
    dsmcCloud& cloud,
    const dictionary& dict
)
:
    dsmcReaction(t, cloud, dict),
    propsDict_(dict.subDict(typeName + "Properties")),
    posMolReactant_(-1),
    productIdsExchange_(),
    exchangeStr_(word::null),
    nTotExchangeReactions_(0),
    nExchangeReactionsPerTimeStep_(0),
    heatOfReactionExchangeJoules_
    (
        readScalar(propsDict_.lookup("heatOfReactionExchange"))
       *physicoChemical::k.value()
    ),
    aCoeff_(readScalar(propsDict_.lookup("aCoeff"))),
    bCoeff_(readScalar(propsDict_.lookup("bCoeff"))),
    volume_(0.0)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

exchangeQK::~exchangeQK()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void exchangeQK::initialConfiguration()
{
    setProperties();

    const word& reactantA = cloud_.typeIdList()[reactantIds_[0]];
    const word& reactantB = cloud_.typeIdList()[reactantIds_[1]];

    const word& productA = cloud_.typeIdList()[productIdsExchange_[0]];
    const word& productB = cloud_.typeIdList()[productIdsExchange_[1]];

    exchangeStr_ = "Exchange reaction " + reactantA + " + "
        + reactantB + " --> " + productA + " + " + productB;
}


bool exchangeQK::tryReactMolecules
(
    const label& typeIdP,
    const label& typeIdQ
)
{
    //- Function used when setting the pair addressing matrix
    const label reactantPId = findIndex(reactantIds_, typeIdP);
    const label reactantQId = findIndex(reactantIds_, typeIdQ);

    //- If both indices were found in the list of reactants, there Ids will be
    //  different from -1
    if ((reactantPId != -1) && (reactantQId != -1))
    {
        //- Case of dissimilar species colliding (by definition)
        if(reactantPId != reactantQId)
        {
            return true;
        }
    }

    return false;
}


void exchangeQK::reaction
(
    dsmcParcel& p,
    dsmcParcel& q,
    const DynamicList<label>& candidateList,
    const List<DynamicList<label> >& candidateSubList,
    const label& candidateP,
    const List<label>& whichSubCell
)
{}


void exchangeQK::reaction(dsmcParcel& p, dsmcParcel& q)
{
    //- Reset the relax switch
    relax_ = true;

    const label typeIdP = p.typeId();
    const label typeIdQ = q.typeId();

    //- Exchange reaction AB + C --> A + BC
    //  If P is the first reactant AB (i.e., not the atom)
    //  NB: Q is necessarily M otherwise this class would not have been selected
    if
    (
        cloud_.constProps(typeIdP).type() != 10
     && cloud_.constProps(typeIdP).type() != 11
    )
    {
        const scalar mP = cloud_.constProps(typeIdP).mass();
        const scalar mQ = cloud_.constProps(typeIdQ).mass();
        const scalar mR = mP*mQ/(mP + mQ);

        const scalar omegaPQ =
            0.5
            *(
                  cloud_.constProps(typeIdP).omega()
                + cloud_.constProps(typeIdQ).omega()
            );

        const scalar cRsqr = magSqr(p.U() - q.U());
        const scalar translationalEnergy = 0.5*mR*cRsqr;

        //- Possible reactions:
        // 1. Exchange reaction

        scalar totalReactionProbability = 0.0;
        scalarList reactionProbabilities(1, 0.0);
        scalarList collisionEnergies(1, 0.0);

        testExchange
        (
            p,
            translationalEnergy,
            omegaPQ,
            collisionEnergies[0],
            totalReactionProbability,
            reactionProbabilities[0]
        );

        //- Decide if an exchange reaction is to occur
        if (totalReactionProbability > cloud_.rndGen().sample01<scalar>())
        {
            exchange(p, q, collisionEnergies[0]);
        }
    }
    else
    {
        //- If P is the second reactant M, then switch arguments in this
        //  function and P will be first
        exchangeQK::reaction(q, p);
    }
}

void exchangeQK::outputResults(const label& counterIndex)
{
    if (writeRatesToTerminal_)
    {
        //- measure density
        const List<DynamicList<dsmcParcel*>>& cellOccupancy = cloud_.cellOccupancy();

        volume_ = 0.0;

        labelList molsReactants(label(2), 0);

        forAll(cellOccupancy, c)
        {
            const List<dsmcParcel*>& parcelsInCell = cellOccupancy[c];

            forAll(parcelsInCell, pIC)
            {
                dsmcParcel* p = parcelsInCell[pIC];

                const label pos = findIndex(reactantIds_, p->typeId());

                if (pos != -1)
                {
                    molsReactants[pos]++;
                }
            }

            volume_ += mesh_.cellVolumes()[c];
        }

        scalar volume = volume_;
        if (Pstream::parRun())
        {
            reduce(volume, sumOp<scalar>());
        }

        scalarList numberDensities(2, cloud_.nParticle()/volume);
        numberDensities[0] *= molsReactants[0];
        numberDensities[1] *= molsReactants[1];

        label nTotExchangeReactions = nTotExchangeReactions_;
        label nExchangeReactionsPerTimeStep = nExchangeReactionsPerTimeStep_;

        const scalar deltaT = mesh_.time().deltaT().value();
        scalar factor = 0.0;

        if (reactantIds_[0] == reactantIds_[1] && numberDensities[0] > 0.0)
        {
            factor = cloud_.nParticle()/
                (
                    counterIndex*deltaT
                   *numberDensities[0]*numberDensities[0]
                   *volume
                );
        }
        else if (numberDensities[0] > 0.0 && numberDensities[1] > 0.0)
        {
            factor = cloud_.nParticle()/
                (
                    counterIndex*deltaT
                   *numberDensities[0]*numberDensities[1]
                   *volume
                );
        }

        if (Pstream::parRun())
        {
            //- Parallel communication
            reduce(molsReactants[0], sumOp<label>());
            reduce(molsReactants[1], sumOp<label>());
            reduce(nTotExchangeReactions, sumOp<label>());
            reduce(nExchangeReactionsPerTimeStep, sumOp<label>());
        }

        const scalar reactionRateExchange = factor*nTotExchangeReactions;

        Info<< exchangeStr_
            << ", reaction rate = " << reactionRateExchange
            << ", nReactions = " << nExchangeReactionsPerTimeStep
            << endl;
    }
    else
    {
        label nTotExchangeReactions = nTotExchangeReactions_;
        label nExchangeReactionsPerTimeStep = nExchangeReactionsPerTimeStep_;

        if (Pstream::parRun())
        {
            //- Parallel communication
            reduce(nTotExchangeReactions, sumOp<label>());
            reduce(nExchangeReactionsPerTimeStep, sumOp<label>());
        }

        if (nTotExchangeReactions > 0)
        {
            Info<< exchangeStr_
                << " is active, nReactions this time step = "
                << nExchangeReactionsPerTimeStep
                << endl;
         }
    }

    nExchangeReactionsPerTimeStep_ = 0;
}

}
// End namespace Foam

// ************************************************************************* //
