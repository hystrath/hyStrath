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
    simplifiedBernoulliTrialsTAS

Description

Stefan Stefanov's simplified Bernoulli trials collision partner selection routine.
Leads to less repeat collisions than Bird's no time counter method.
See SIAM J. Sci. Comput. 33(2) 667-702 for details.
We do not use the half timestep implementation.

\*----------------------------------------------------------------------------*/

#include "simplifiedBernoulliTrialsTAS.H"
#include "addToRunTimeSelectionTable.H"
#include "boundBox.H"

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(simplifiedBernoulliTrialsTAS, 0);

addToRunTimeSelectionTable(collisionPartnerSelection, simplifiedBernoulliTrialsTAS, dictionary);



// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //
pointField simplifiedBernoulliTrialsTAS::cellPoints(const label& cell)
{
    const labelList& points = mesh_.cellPoints()[cell];

    pointField vectorPoints(points.size(), vector::zero);

    forAll(points, p)
    {
        vectorPoints[p] = mesh_.points()[points[p]];
    }

    return vectorPoints;
}

void simplifiedBernoulliTrialsTAS::checkSubCelling()
{

}

void simplifiedBernoulliTrialsTAS::measureLocalDensity()
{
    const List< DynamicList<dsmcParcel*> >& cellOccupancy
        = cloud_.cellOccupancy();
    
    forAll(cellOccupancy, cell)
    {
        const List<dsmcParcel*>& parcelsInCell = cellOccupancy[cell];
        rhoN_[cell] = parcelsInCell.size()/mesh_.cellVolumes()[cell];
    }

}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

//- Construct from components
simplifiedBernoulliTrialsTAS::simplifiedBernoulliTrialsTAS
(
    const polyMesh& mesh,
    dsmcCloud& cloud,
    const dictionary& dict
)
:
    collisionPartnerSelection(mesh, cloud, dict),
    propsDict_(dict.subDict(typeName + "Properties")),
    n_(readLabel(propsDict_.lookup("nParcelsPerSubCell"))),
    nSteps_(readLabel(propsDict_.lookup("nTimeSteps"))),
    counter_(0),
    rhoN_(mesh.nCells(), 0.0),
    startPoints_(mesh.nCells(), vector::zero),
    nSlices_(mesh.nCells()),
    nSubCells_(mesh.nCells(), -1),
    binWidths_(mesh.nCells(), vector::zero)
{
    forAll(nSlices_, c)
    {
        nSlices_[c].setSize(3, 0);
    }

    if(n_ <= 0)
    {
        FatalErrorIn("simplifiedBernoulliTrialsTAS::simplifiedBernoulliTrialsTAS") << nl
            << "nParcelsPerSubCell = " << n_
            << nl << "ERROR: Choose a number greater than 0."
            << nl << abort(FatalError);
    }
}



// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

simplifiedBernoulliTrialsTAS::~simplifiedBernoulliTrialsTAS()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void simplifiedBernoulliTrialsTAS::initialConfiguration()
{
    measureLocalDensity();

    // 1. Create boundbox around each cell and find the starting point
    //    i.e. the minimum point at the base of the bound box

    // 2. Identify the number of slices in the x,y z coordinates and their binWidths

    for (label c = 0; c < mesh_.nCells(); c++)
    {
        boundBox bb(cellPoints(c), false);

        startPoints_[c] = bb.min();

        scalar bbVol = bb.span().x() * bb.span().y() * bb.span().z();
        scalar nTot = rhoN_[c]*bbVol/n_;

        // test for no subcelling
        if(nTot <= 1.0)
        {
            nSlices_[c][0] = 1;
            nSlices_[c][1] = 1;
            nSlices_[c][2] = 1;
        }
        else
        {
            scalar deltaLCubic = Foam::pow(bbVol, (1.0/3.0) );
            scalar spacingCubic = deltaLCubic/Foam::pow(nTot, (1.0/3.0) );
    
            scalar nX = bb.span().x()/spacingCubic;
            scalar nY = bb.span().y()/spacingCubic;
            scalar nZ = bb.span().z()/spacingCubic;
    
            nSlices_[c][0] = label(nX+0.5);
            nSlices_[c][1] = label(nY+0.5);
            nSlices_[c][2] = label(nZ+0.5);
    
            // test for zero slices (1 is the bare minimum) and modfiy
            if(nSlices_[c][0] == 0)
            {
                nSlices_[c][0] = 1;
            }
            if(nSlices_[c][1] == 0)
            {
                nSlices_[c][1] = 1;
            }
            if(nSlices_[c][2] == 0)
            {
                nSlices_[c][2] = 1;
            }
        }

        binWidths_[c].x() = bb.span().x()/scalar(nSlices_[c][0]);
        binWidths_[c].y() = bb.span().y()/scalar(nSlices_[c][1]);
        binWidths_[c].z() = bb.span().z()/scalar(nSlices_[c][2]);

        nSubCells_[c] = nSlices_[c][0]*nSlices_[c][1]*nSlices_[c][2];
    }
}


void simplifiedBernoulliTrialsTAS::collide()
{
	
	// Measure instantaneous density and reset sub-cells
    if(counter_ >= nSteps_)
    {
        initialConfiguration();
        counter_ = 0;
    }

    counter_++;

	
	if (!cloud_.binaryCollision().active())
    {
        return;
    }

    scalar deltaT = cloud_.mesh().time().deltaTValue();

    label collisionCandidates = 0;

    label collisions = 0;
	
	const List<DynamicList<dsmcParcel*> > cellOccupancy = cloud_.cellOccupancy();

    const polyMesh& mesh = cloud_.mesh();
	label nX;
    label nY;
    label nZ;
    label subCell;
    vector pS;

    forAll(cellOccupancy, cellI)
    {
        const DynamicList<dsmcParcel*>& cellParcels(cellOccupancy[cellI]);

        label nC(cellParcels.size());

        if (nC > 1)
        {

			const point& cC = mesh.cellCentres()[cellI];

            List<DynamicList<label> > subCells(nSubCells_[cellI]);

            forAll(cellParcels, i)
            {
                const dsmcParcel& p = *cellParcels[i];
                pS = p.position() - startPoints_[cellI];
                nX = label((pS & vector(1, 0, 0))/binWidths_[cellI].x());
                nY = label((pS & vector(0, 1, 0))/binWidths_[cellI].y());
                nZ = label((pS & vector(0, 0, 1))/binWidths_[cellI].z());
        
                subCell = nX + nY*nSlices_[cellI][0] + nZ*nSlices_[cellI][0]*nSlices_[cellI][1];
        
                subCells[subCell].append(i);
        
            }

        
            scalar prob1 = (cloud_.nParticle()*deltaT)/(mesh.cellVolumes()[cellI]/nSubCells_[cellI]);
            label k = -1;
            label candidateP = -1;
            label candidateQ = -1;
        
        // loop over sub cells
            forAll(subCells, i)
            {
                subCells[i].shrink();
                label nCS = subCells[i].size();
                
                for(label p = 0 ; p < nCS-1 ; p++)
                {                
                    // Select the first collision candidate
                    candidateP = p;
                
                    k = nCS - p - 1;
                    candidateQ = p + rndGen_.integer(1, k);
                    
                    dsmcParcel& parcelP = *cellParcels[subCells[i][candidateP]];
                    dsmcParcel& parcelQ = *cellParcels[subCells[i][candidateQ]];

                    scalar sigmaTcR = cloud_.binaryCollision().sigmaTcR
                    (
                        parcelP,
                        parcelQ
                    );
                
                    scalar Probability = k*prob1*sigmaTcR*0.5;

    // 					if (Probability > rndGen_.scalar01())
    // 					{
    // 						binaryCollision().collide
    // 						(
    // 							parcelP,
    // 							parcelQ
    // 						);
    // 
    // 						collisions++;
    // 					}
                    if (Probability > rndGen_.scalar01())
                    {
                        // chemical reactions

                        // find which reaction model parcel p and q should use
                        label rMId = cloud_.reactions().returnModelId(parcelP, parcelQ);

    //                             Info << " parcelP id: " <<  parcelP.typeId() 
    //                                 << " parcelQ id: " << parcelQ.typeId()
    //                                 << " reaction model: " << rMId
    //                                 << endl;

                        if(rMId != -1)
                        {
                            // try to react molecules
                            if(cloud_.reactions().reactions()[rMId]->reactWithLists())
                            {
                                // so far for recombination only
    //                                     reactions_.reactions()[rMId]->reaction
    //                                     (
    //                                         parcelP,
    //                                         parcelQ,
    //                                         candidateList,
    //                                         candidateSubList,
    //                                         candidateP,
    //                                         whichSubCell
    //                                     );
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

    //                                 buildCellOccupancy();
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
