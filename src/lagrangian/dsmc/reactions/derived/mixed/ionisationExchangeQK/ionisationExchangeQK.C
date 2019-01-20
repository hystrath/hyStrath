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

#include "ionisationExchangeQK.H"
#include "addToRunTimeSelectionTable.H"
#include "fvc.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(ionisationExchangeQK, 0);

addToRunTimeSelectionTable
(
    dsmcReaction,
    ionisationExchangeQK,
    dictionary
);


// * * * * * * * * * * *  Protected Member functions * * * * * * * * * * * * //

void ionisationExchangeQK::setProperties()
{
    ionisationQK::setProperties();
    exchangeQK::setProperties();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
ionisationExchangeQK::ionisationExchangeQK
(
    Time& t,
    dsmcCloud& cloud,
    const dictionary& dict
)
:
    dsmcReaction(t, cloud, dict),
    ionisationQK(t, cloud, dict),
    exchangeQK(t, cloud, dict)
    //propsDict_(dict.subDict(typeName + "Properties")),
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

ionisationExchangeQK::~ionisationExchangeQK()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void ionisationExchangeQK::initialConfiguration()
{
    ionisationQK::initialConfiguration();
    exchangeQK::initialConfiguration();
}


bool ionisationExchangeQK::tryReactMolecules
(
    const label& typeIdP,
    const label& typeIdQ
)
{
    return ionisationQK::tryReactMolecules(typeIdP, typeIdQ);
}


void ionisationExchangeQK::reaction
(
    dsmcParcel& p,
    dsmcParcel& q,
    const DynamicList<label>& candidateList,
    const List<DynamicList<label> >& candidateSubList,
    const label& candidateP,
    const List<label>& whichSubCell
)
{}


void ionisationExchangeQK::reaction(dsmcParcel& p, dsmcParcel& q)
{
    //- Reset the relax switch
    relax_ = true;
    
    const label typeIdP = p.typeId();
    const label typeIdQ = q.typeId();
    
    if (typeIdP == reactantIds_[0]) 
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
        // 1. Ionisation of P
        // 2. Ionisation of Q
        // 3. Exchange
        
        scalar totalReactionProbability = 0.0;
        scalarList reactionProbabilities(3, 0.0);
        scalarList collisionEnergies(3, 0.0);
        
        ionisationQK::testIonisation
        (
            p,
            translationalEnergy,
            0,
            collisionEnergies[0],
            totalReactionProbability,
            reactionProbabilities[0]
        );
        
        ionisationQK::testIonisation
        (
            q,
            translationalEnergy,
            1,
            collisionEnergies[1],
            totalReactionProbability,
            reactionProbabilities[1]
        );
        
        if (exchangeQK::posMolReactant_ == 0)
        {
            exchangeQK::testExchange
            (
                p,
                translationalEnergy,
                omegaPQ,
                collisionEnergies[2],
                totalReactionProbability,
                reactionProbabilities[2]
            );
        }
        else
        {
            exchangeQK::testExchange
            (
                q,
                translationalEnergy,
                omegaPQ,
                collisionEnergies[2],
                totalReactionProbability,
                reactionProbabilities[2]
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
                            ionisationQK::ioniseParticleByPartner
                            (
                                p, q, 0, collisionEnergies[i]
                            );
                            //- There can't be another reaction: break
                            break;
                        }
                        
                        if (i == 1)
                        {
                            //- Ionisation of Q is to occur
                            ionisationQK::ioniseParticleByPartner
                            (
                                q, p, 1, collisionEnergies[i]
                            );
                            //- There can't be another reaction: break
                            break;
                        }
                        
                        if (i == 2)
                        {
                            //- Exchange reaction
                            if (exchangeQK::posMolReactant_ == 0)
                            {
                                exchangeQK::exchange
                                (
                                    p, q, collisionEnergies[i]
                                );
                            }
                            else
                            {
                                exchangeQK::exchange
                                (
                                    q, p, collisionEnergies[i]
                                );
                            }
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
        //  If P is the second reactant, then switch arguments in this
        //  function and P will be first
        ionisationExchangeQK::reaction(q, p);
    }
}


inline label ionisationExchangeQK::nReactionsPerTimeStep() const
{
    return ionisationQK::nReactionsPerTimeStep()
        + exchangeQK::nReactionsPerTimeStep();
}


void ionisationExchangeQK::outputResults
(
    const label& counterIndex
)
{
    if (writeRatesToTerminal_)
    {
        ionisationQK::outputResults(counterIndex);
        exchangeQK::outputResults(counterIndex);
    }
}

}
// End namespace Foam

// ************************************************************************* //
