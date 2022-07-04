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

#include "ionisationQK.H"
#include "addToRunTimeSelectionTable.H"
#include "fvc.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(ionisationQK, 0);

addToRunTimeSelectionTable(dsmcReaction, ionisationQK, dictionary);


// * * * * * * * * * * *  Protected Member functions * * * * * * * * * * * * //

void ionisationQK::setProperties()
{
    dsmcReaction::setProperties();

    if (reactantIds_.size() != 2)
    {
        //- There must be exactly 2 reactants
        FatalErrorIn("ionisationQK::setProperties()")
            << "For reaction named " << reactionName_ << nl
            << "There should be two reactants, instead of "
            << reactantIds_.size() << nl
            << exit(FatalError);
    }

    ionisationStr_.setSize(reactantIds_.size(), word::null);
    nTotIonisationReactions_.setSize(reactantIds_.size(), 0);
    nIonisationReactionsPerTimeStep_.setSize(reactantIds_.size(), 0);
    heatOfReactionIonisationJoules_.setSize(reactantIds_.size(), 0.0);

    forAll(reactantIds_, r)
    {
        //- Check that this reactant is not an electron
        if (reactantTypes_[r] != 0)
        {
            heatOfReactionIonisationJoules_[r] = physicoChemical::k.value()
                *cloud_.constProps(reactantIds_[r]).ionisationTemperature();
        }

        if (reactantTypes_[r] >= 30)
        {
            //- Polyatomic molecules are not handled yet
            FatalErrorIn("ionisationQK::setProperties()")
                << "For reaction named " << reactionName_ << nl
                << "Reactions are currently only implemented for monatomic and diatomic species"
                << " This is a polyatomic:"
                << cloud_.typeIdList()[reactantIds_[r]] << nl
                << exit(FatalError);
        }
    }

    //- Reading in ionisation products
    const List<wordList> productsIonisation(propsDict_.lookup("ionisationProducts"));

    if (productsIonisation.size() != 2)
    {
        FatalErrorIn("ionisationQK::setProperties()")
            << "For reaction named " << reactionName_ << nl
            << "There should be two lists of products, instead of "
            << productsIonisation.size() << nl
            << "NB: a list can be left empty" << nl
            << exit(FatalError);
    }

    productIdsIonisation_.setSize(productsIonisation.size());

    forAll(productIdsIonisation_, r)
    {
        if (productsIonisation[r].size() > 0)
        {
            //- Check that there is no ionisation products should the reactant
            //  be an electron
            if (reactantTypes_[r] == 0)
            {
                FatalErrorIn("ionisationQK::setProperties()")
                    << "For reaction named " << reactionName_ << nl
                    << "Reactant " << cloud_.typeIdList()[reactantIds_[r]]
                    << " is not a molecule " << nl
                    << "and therefore, there should be no ionisation products "
                    << "instead of " << productsIonisation[r]
                    << exit(FatalError);
            }

            //- Check that there are two products
            if (productsIonisation[r].size() != 2)
            {
                FatalErrorIn("ionisationQK::setProperties()")
                    << "For reaction named " << reactionName_ << nl
                    << "There should be 2 ionisation products for molecule "
                    << cloud_.typeIdList()[reactantIds_[r]] << " instead of "
                    << productsIonisation[r].size() << ", that is "
                    << productsIonisation[r]
                    << exit(FatalError);
            }
        }
        else
        {
            //- If this species is not charged, it should have ionisation
            //  products
            if (cloud_.constProps(reactantIds_[r]).charge() == 0)
            {
                FatalErrorIn("ionisationQK::setProperties()")
                    << "For reaction named " << reactionName_ << nl
                    << "Reactant " << cloud_.typeIdList()[reactantIds_[r]]
                    << " is a molecule " << nl
                    << "and therefore, it should have ionisation products "
                    << "instead of " << productsIonisation[r]
                    << exit(FatalError);
            }
        }

        productIdsIonisation_[r].setSize(productsIonisation[r].size());

        forAll(productIdsIonisation_[r], p)
        {
            productIdsIonisation_[r][p] =
                findIndex
                (
                    cloud_.typeIdList(),
                    productsIonisation[r][p]
                );

            //- Check that products belong to the typeIdList as defined in
            //  constant/dsmcProperties
            if (productIdsIonisation_[r][p] == -1)
            {
                FatalErrorIn("ionisationQK::setProperties()")
                    << "For reaction named " << reactionName_ << nl
                    << "Cannot find type id: " << productsIonisation[r][p] << nl
                    << exit(FatalError);
            }

            //- Check that product is a charged particle
            if (cloud_.constProps(productIdsIonisation_[r][p]).charge() == 0)
            {
                FatalErrorIn("ionisationQK::setProperties()")
                    << "For reaction named " << reactionName_ << nl
                    << "Ionisation products must be charged particles: "
                    << productsIonisation[r][p] << nl
                    << exit(FatalError);
            }
        }
    }
}


void ionisationQK::testIonisation
(
    const dsmcParcel& p,
    const scalar translationalEnergy,
    const label nReac,
    scalar& collisionEnergy,
    scalar& totalReactionProbability,
    scalar& reactionProbability
)
{
    const label typeIdP = p.typeId();

    if (cloud_.constProps(typeIdP).type() != 0)
    {
        const scalar EEleP =
            cloud_.constProps(typeIdP).electronicEnergyList()[p.ELevel()];
        collisionEnergy = translationalEnergy + EEleP;

        //- Condition for the ionisation of P
        if (collisionEnergy > heatOfReactionIonisationJoules_[nReac])
        {
            //- Add reaction to the list of competing reactions with
            //  probability reactionProbability
            reactionProbability = 1.0;
            totalReactionProbability += reactionProbability;
        }
    }
}


void ionisationQK::ioniseParticleByPartner
(
    dsmcParcel& p,
    dsmcParcel& q,
    const label nR,
    scalar collisionEnergy
)
{
    const label typeIdP = p.typeId();
    const label typeIdQ = q.typeId();

    label nReac = nR;
    if (typeIdP == typeIdQ)
    {
        nReac = 0;
    }

    //- Ionisation of parcel p
    nTotIonisationReactions_[nReac]++;
    nIonisationReactionsPerTimeStep_[nReac]++;

    if (allowSplitting_)
    {
        relax_ = false;

        //- The collision energy is being subtracted the heat of reaction
        collisionEnergy -= heatOfReactionIonisationJoules_[nR];

        const scalar omegaPQ =
            0.5
            *(
                  cloud_.constProps(typeIdP).omega()
                + cloud_.constProps(typeIdQ).omega()
            );

        //- Energy redistribution for particle Q
        cloud_.binaryCollision().redistribute
        (
            q, collisionEnergy, omegaPQ, true
        );

        const scalar mP = cloud_.constProps(typeIdP).mass();
        const scalar mQ = cloud_.constProps(typeIdQ).mass();
        const scalar mR = mP*mQ/(mP + mQ);

        scalar relVelNonIonisedParticle = sqrt(2.0*collisionEnergy/mR);

        //- Post-collision velocities for P (pre-reaction) and Q
        cloud_.binaryCollision().postCollisionVelocities
        (
            typeIdP,
            typeIdQ,
            p.U(),
            q.U(),
            relVelNonIonisedParticle
        );

        //- Post-reaction
        const label typeId1 = productIdsIonisation_[nR][0];
        const label typeId2 = productIdsIonisation_[nR][1];

        //- Mass of products 1 and 2
        const scalar mP1 = cloud_.constProps(typeId1).mass();
        const scalar mP2 = cloud_.constProps(typeId2).mass();
        const scalar mRproducts = mP1*mP2/(mP1 + mP2);

        //- Energy left for the 2 products
        const scalar ERotP = p.ERot();
        const scalar EVibP = cloud_.constProps(typeIdP).eVib_tot(p.vibLevel());
        //- Assumption: no energy redistribution for the particle being split
        //  All the remaining energy is stored in the translational mode
        const scalar translationalEnergyLeft = ERotP + EVibP;
        const scalar cRproducts = sqrt(2.0*translationalEnergyLeft/mRproducts);

        vector UP2 = vector::zero;
        cloud_.binaryCollision().postReactionVelocities
        (
            typeId1,
            typeId2,
            p.U(),
            UP2,
            cRproducts
        );

        // Molecule P ionises
        // The molecule is first transformed into the ionised particle ...
        p.typeId() = typeId1;
        p.vibLevel() = 0;
        p.ERot() = 0.0;
        p.ELevel() = 0;

        // ... and an electron is then inserted in the same cell
        const scalar product2_ERot = 0.0;
        const labelList product2_vibLevel
        (
            cloud_.constProps(typeId2).nVibrationalModes(),
            0
        );
        const label product2_ELevel = 0;
        cloud_.addNewParcel
        (
            p.position(),
            UP2,
            p.RWF(),
            product2_ERot,
            product2_ELevel,
            p.cell(),
            p.tetFace(),
            p.tetPt(),
            typeId2,
            -1,
            p.classification(),
            product2_vibLevel
        );
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
ionisationQK::ionisationQK
(
    Time& t,
    dsmcCloud& cloud,
    const dictionary& dict
)
:
    dsmcReaction(t, cloud, dict),
    propsDict_(dict.subDict(typeName + "Properties")),
    productIdsIonisation_(),
    ionisationStr_(),
    nTotIonisationReactions_(),
    nIonisationReactionsPerTimeStep_(),
    heatOfReactionIonisationJoules_(),
    volume_(0.0)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

ionisationQK::~ionisationQK()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void ionisationQK::initialConfiguration()
{
    setProperties();

    const word& reactantA = cloud_.typeIdList()[reactantIds_[0]];
    const word& reactantB = cloud_.typeIdList()[reactantIds_[1]];

    if (productIdsIonisation_[0].size())
    {
        const word& productA = cloud_.typeIdList()[productIdsIonisation_[0][0]];
        const word& productB = cloud_.typeIdList()[productIdsIonisation_[0][1]];
        ionisationStr_[0] = "Ionisation reaction " + reactantA + " + "
            + reactantB + " --> " + productA + " + " + productB + " + "
            + reactantB;
    }
    if (productIdsIonisation_[1].size() and reactantA != reactantB)
    {
        const word& productC = cloud_.typeIdList()[productIdsIonisation_[1][0]];
        const word& productD = cloud_.typeIdList()[productIdsIonisation_[1][1]];

        ionisationStr_[1] = "Ionisation reaction " + reactantB + " + "
            + reactantA + " --> " + productC + " + " + productD + " + "
            + reactantA;
    }
}


bool ionisationQK::tryReactMolecules
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
        //- Case of similar species colliding
        if((reactantPId == reactantQId) and (reactantIds_[0] == reactantIds_[1]))
        {
            return true;
        }

        //- Case of dissimilar species colliding
        if((reactantPId != reactantQId) and (reactantIds_[0] != reactantIds_[1]))
        {
            return true;
        }
    }

    return false;
}


void ionisationQK::reaction
(
    dsmcParcel& p,
    dsmcParcel& q,
    const DynamicList<label>& candidateList,
    const List<DynamicList<label> >& candidateSubList,
    const label& candidateP,
    const List<label>& whichSubCell
)
{}


void ionisationQK::reaction(dsmcParcel& p, dsmcParcel& q)
{
    //- Reset the relax switch
    relax_ = true;

    const label typeIdP = p.typeId();
    const label typeIdQ = q.typeId();

    //- Ionisation reaction A + M --> A+ + E- + M
    //  If P is the first reactant A
    //  NB: Q is necessarily M otherwise this class would not have been selected
    if (typeIdP == reactantIds_[0])
    {
        const scalar mP = cloud_.constProps(typeIdP).mass();
        const scalar mQ = cloud_.constProps(typeIdQ).mass();
        const scalar mR = mP*mQ/(mP + mQ);

        const scalar cRsqr = magSqr(p.U() - q.U());
        const scalar translationalEnergy = 0.5*mR*cRsqr;

        //- Possible reactions:
        // 1. Ionisation of P
        // 2. Ionisation of Q

        scalar totalReactionProbability = 0.0;
        scalarList reactionProbabilities(2, 0.0);
        scalarList collisionEnergies(2, 0.0);

        testIonisation
        (
            p,
            translationalEnergy,
            0,
            collisionEnergies[0],
            totalReactionProbability,
            reactionProbabilities[0]
        );

        if (typeIdP != typeIdQ)
        {
            testIonisation
            (
                q,
                translationalEnergy,
                1,
                collisionEnergies[1],
                totalReactionProbability,
                reactionProbabilities[1]
            );
        }

        //- Decide if a reaction is to occur
        if (totalReactionProbability > cloud_.rndGen().sample01<scalar>())
        {
            //- A chemical reaction is to occur, normalise probabilities
            const scalarList normalisedProbabilities =
                reactionProbabilities/totalReactionProbability;

            //- Sort normalised probability indices in decreasing order
            //  for identical probabilities, random shuffle
            const labelList sortedNormalisedProbabilityIndices =
                decreasing_sort_indices(normalisedProbabilities);
            scalar cumulativeProbability = 0.0;

            forAll(sortedNormalisedProbabilityIndices, idx)
            {
                const label i = sortedNormalisedProbabilityIndices[idx];

                //- If current reaction can't occur, end the search
                if (normalisedProbabilities[i] > SMALL)
                {
                    cumulativeProbability += normalisedProbabilities[i];

                    if (cumulativeProbability > cloud_.rndGen().sample01<scalar>())
                    {
                        //- Current reaction is to occur
                        if (i == 0)
                        {
                            //- Ionisation of P is to occur
                            ioniseParticleByPartner
                            (
                                p, q, i, collisionEnergies[i]
                            );
                            //- There can't be another reaction: break
                            break;
                        }

                        if (i == 1)
                        {
                            //- Ionisation of Q is to occur
                            ioniseParticleByPartner
                            (
                                q, p, i, collisionEnergies[i]
                            );
                            //- There can't be another reaction: break
                            break;
                        }
                    }
                }
                else
                {
                    //- All the following possible reactions have a probability
                    //  of zero
                    break;
                }
            }
        }
    }
    else
    {
        //- Ionisation reaction A + M --> A+ + E- + M
        //  If P is the second reactant M, then switch arguments in this
        //  function and P will be first
        ionisationQK::reaction(q, p);
    }
}

void ionisationQK::outputResults(const label& counterIndex)
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

        labelList nTotIonisationReactions = nTotIonisationReactions_;
        labelList nIonisationReactionsPerTimeStep = nIonisationReactionsPerTimeStep_;

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

        for (label k=0; k<2; k++)
        {
            if (ionisationStr_[k].size())
            {
                if (Pstream::parRun())
                {
                    //- Parallel communication
                    reduce(molsReactants[k], sumOp<label>());
                    reduce(nTotIonisationReactions[k], sumOp<label>());
                    reduce(nIonisationReactionsPerTimeStep[k], sumOp<label>());
                }

                const scalar reactionRateIonisation = factor*nTotIonisationReactions[k];

                Info<< ionisationStr_[k]
                    << ", reaction rate = " << reactionRateIonisation
                    << ", nReactions = " << nIonisationReactionsPerTimeStep[k]
                    << endl;
            }
        }
    }
    else
    {
        labelList nTotIonisationReactions = nTotIonisationReactions_;
        labelList nIonisationReactionsPerTimeStep = nIonisationReactionsPerTimeStep_;

        for (label k=0; k<2; k++)
        {
            if (ionisationStr_[k].size())
            {
                if (Pstream::parRun())
                {
                    //- Parallel communication
                    reduce(nTotIonisationReactions[k], sumOp<label>());
                    reduce(nIonisationReactionsPerTimeStep[k], sumOp<label>());
                }

                if (nTotIonisationReactions[k] > 0)
                {
                    Info<< ionisationStr_[k]
                        << " is active, nReactions this time step = "
                        << nIonisationReactionsPerTimeStep[k]
                        << endl;
                 }
             }
        }
    }

    nIonisationReactionsPerTimeStep_ = 0;
}

}
// End namespace Foam

// ************************************************************************* //
