/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2005 OpenCFD Ltd.
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

Class
    simplifiedBernoulliTrials

Description

Stefan Stefanov's simplified Bernoulli trials collision partner selection routine.
Leads to less repeat collisions than Bird's no time counter method.
See SIAM J. Sci. Comput. 33(2) 667-702 for details.
We do not use the half timestep implementation.

\*----------------------------------------------------------------------------*/

#include "simplifiedBernoulliTrials.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(simplifiedBernoulliTrials, 0);

addToRunTimeSelectionTable
(
    collisionPartnerSelection,
    simplifiedBernoulliTrials,
    dictionary
);

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

//- Construct from components
simplifiedBernoulliTrials::simplifiedBernoulliTrials
(
    const polyMesh& mesh,
    dsmcCloud& cloud,
    const dictionary& dict
)
:
    collisionPartnerSelection(mesh, cloud, dict)
{}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

simplifiedBernoulliTrials::~simplifiedBernoulliTrials()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void simplifiedBernoulliTrials::initialConfiguration()
{

}


void simplifiedBernoulliTrials::collide()
{
    if (!cloud_.binaryCollision().active())
    {
        return;
    }

    //- Temporary storage for subCells
    List<DynamicList<label> > subCells(8);

    label collisionCandidates = 0;

    label collisions = 0;
	
    const List<DynamicList<dsmcParcel*> > cellOccupancy = cloud_.cellOccupancy();

    const polyMesh& mesh = cloud_.mesh();

    forAll(cellOccupancy, cellI)
    {
        const scalar deltaT = cloud_.deltaTValue(cellI);
        
        const scalar nParticle = cloud_.nParticles(cellI, true);
        
        const DynamicList<dsmcParcel*>& cellParcels(cellOccupancy[cellI]);

        label nC(cellParcels.size());

        if (nC > 1)
        {   
            scalar prob1 = (nParticle*deltaT)/(mesh.cellVolumes()[cellI]);
            label k = -1;
            label candidateP = -1;
            label candidateQ = -1;
                
            for(label p = 0 ; p < nC-1 ; p++)
            {                
                // Select the first collision candidate
                candidateP = p;
            
                k = nC-1 - p;
                //label random = rndGen_.position<label>(1, k); OLD
                label random = 1 + rndGen_.sample01<scalar>()*(k+1);
                candidateQ = p + random;
                
                dsmcParcel& parcelP = *cellParcels[candidateP];
                dsmcParcel& parcelQ = *cellParcels[candidateQ];  

                scalar sigmaTcR = cloud_.binaryCollision().sigmaTcR
                (
                    parcelP,
                    parcelQ
                );
            
                scalar Probability = k*prob1*sigmaTcR;

                if (Probability > rndGen_.sample01<scalar>())
                {
                    // chemical reactions

                    // find which reaction model parcel p and q should use
                    label rMId = cloud_.reactions().returnModelId(parcelP, parcelQ);

                    if(rMId != -1)
                    {
                        // try to react molecules
                        if(cloud_.reactions().reactions()[rMId]->reactWithLists())
                        {

                        }
                        else
                        {
                            cloud_.reactions().reactions()[rMId]->reaction
                            (
                                parcelP,
                                parcelQ
                            );                                    
                        }
                        // if reaction unsuccessful use conventional collision model
                        if(cloud_.reactions().reactions()[rMId]->relax())
                        {
                            cloud_.binaryCollision().collide
                            (
                                parcelP,
                                parcelQ,
                                cellI
                            );
                        }
                    }
                    else // if reaction model not found, use conventional collision model
                    {
                        cloud_.binaryCollision().collide
                        (
                            parcelP,
                            parcelQ,
                            cellI
                        );
                    }
                    
                    collisions++;
                }
            }
        }
    }

    reduce(collisions, sumOp<label>());

    reduce(collisionCandidates, sumOp<label>());

//     sigmaTcRMax_.correctBoundaryConditions();

    if (collisions > 0)
    {
        Info<< "    Collisions                      = "
            << collisions << endl;
    }
    else
    {
        Info<< "    No collisions" << endl;
    }
}

// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //



// * * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
