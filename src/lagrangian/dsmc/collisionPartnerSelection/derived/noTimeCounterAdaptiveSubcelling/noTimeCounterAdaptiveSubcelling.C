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
    noTimeCounterAdaptiveSubcelling

Description

\*----------------------------------------------------------------------------*/

#include "noTimeCounterAdaptiveSubcelling.H"
#include "addToRunTimeSelectionTable.H"
#include "boundBox.H"

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(noTimeCounterAdaptiveSubcelling, 0);

addToRunTimeSelectionTable(collisionPartnerSelection, noTimeCounterAdaptiveSubcelling, dictionary);



// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

pointField noTimeCounterAdaptiveSubcelling::cellPoints(const label& cell)
{
    const labelList& points = mesh_.cellPoints()[cell];

    pointField vectorPoints(points.size(), vector::zero);

    forAll(points, p)
    {
        vectorPoints[p] = mesh_.points()[points[p]];
    }

    return vectorPoints;
}

void noTimeCounterAdaptiveSubcelling::checkSubCelling()
{

}

void noTimeCounterAdaptiveSubcelling::measureLocalDensity()
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
noTimeCounterAdaptiveSubcelling::noTimeCounterAdaptiveSubcelling
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
        FatalErrorIn("noTimeCounterAdaptiveSubcelling::noTimeCounterAdaptiveSubcelling") << nl
            << "nParcelsPerSubCell = " << n_
            << nl << "ERROR: Choose a number greater than 0."
            << nl << abort(FatalError);
    }
}



// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

noTimeCounterAdaptiveSubcelling::~noTimeCounterAdaptiveSubcelling()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void noTimeCounterAdaptiveSubcelling::initialConfiguration()
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

void noTimeCounterAdaptiveSubcelling::collide()
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
        const scalar deltaT = cloud_.deltaTValue(cellI);
        
        const DynamicList<dsmcParcel*>& cellParcels(cellOccupancy[cellI]);

        label nC(cellParcels.size());

        if (nC > 1)
        {
            // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

            List<DynamicList<label> > subCells(nSubCells_[cellI]);

            List<label> whichSubCell(cellParcels.size());
        
            forAll(cellParcels, i)
            {
                const dsmcParcel& p = *cellParcels[i];
                pS = p.position() - startPoints_[cellI];
                nX = label((pS & vector(1, 0, 0))/binWidths_[cellI].x());
                nY = label((pS & vector(0, 1, 0))/binWidths_[cellI].y());
                nZ = label((pS & vector(0, 0, 1))/binWidths_[cellI].z());
        
                subCell = nX + nY*nSlices_[cellI][0] + nZ*nSlices_[cellI][0]*nSlices_[cellI][1];
        
                subCells[subCell].append(i);
        
                whichSubCell[i] = subCell;
            }

            // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

            scalar sigmaTcRMax = cloud_.sigmaTcRMax()[cellI];

            scalar selectedPairs =
                cloud_.collisionSelectionRemainder()[cellI]
                + 0.5*nC*(nC - 1)*cloud_.nParticles(cellI, true)*sigmaTcRMax*deltaT
                /mesh.cellVolumes()[cellI];
               
            label nCandidates(selectedPairs);

            cloud_.collisionSelectionRemainder()[cellI] = selectedPairs - nCandidates;

            collisionCandidates += nCandidates;

            for (label c = 0; c < nCandidates; c++)
            {
                // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                // subCell candidate selection procedure

                // Select the first collision candidate
                label candidateP = rndGen_.integer(0, nC - 1);

                // Declare the second collision candidate
                label candidateQ = -1;

                const List<label>& subCellPs = subCells[whichSubCell[candidateP]];

                label nSC = subCellPs.size();

                if (nSC > 1)
                {
                    // If there are two or more particle in a subCell, choose
                    // another from the same cell.  If the same candidate is
                    // chosen, choose again.

                    do
                    {
                        candidateQ = subCellPs[rndGen_.integer(0, nSC - 1)];

                    } while (candidateP == candidateQ);
                }
                else
                {
                    // Select a possible second collision candidate from the
                    // whole cell.  If the same candidate is chosen, choose
                    // again.

                    do
                    {
                        candidateQ = rndGen_.integer(0, nC - 1);

                    } while (candidateP == candidateQ);
                }

                // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                // uniform candidate selection procedure

                // // Select the first collision candidate
                // label candidateP = rndGen_.integer(0, nC-1);

                // // Select a possible second collision candidate
                // label candidateQ = rndGen_.integer(0, nC-1);

                // // If the same candidate is chosen, choose again
                // while (candidateP == candidateQ)
                // {
                //     candidateQ = rndGen_.integer(0, nC-1);
                // }

                // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

                dsmcParcel& parcelP = *cellParcels[candidateP];
                dsmcParcel& parcelQ = *cellParcels[candidateQ];

                scalar sigmaTcR = cloud_.binaryCollision().sigmaTcR
                (
                    parcelP,
                    parcelQ
                );

                // Update the maximum value of sigmaTcR stored, but use the
                // initial value in the acceptance-rejection criteria because
                // the number of collision candidates selected was based on this

                if (sigmaTcR > cloud_.sigmaTcRMax()[cellI])
                {
                    cloud_.sigmaTcRMax()[cellI] = sigmaTcR;
                }

                if ((sigmaTcR/sigmaTcRMax) > rndGen_.scalar01())
                {
                    // chemical reactions

                    // find which reaction model parcel p and q should use
                    label rMId = cloud_.reactions().returnModelId(parcelP, parcelQ);

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

    cloud_.sigmaTcRMax().correctBoundaryConditions();

    if (collisionCandidates)
    {
        Info<< "    Collisions                      = "
            << collisions << nl
            << "    Acceptance rate                 = "
            << scalar(collisions)/scalar(collisionCandidates) << nl
            << endl;
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
