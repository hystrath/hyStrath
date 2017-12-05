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
    noTimeCounterSubCycled

Description

\*----------------------------------------------------------------------------*/

#include "noTimeCounterSubCycled.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(noTimeCounterSubCycled, 0);

addToRunTimeSelectionTable(collisionPartnerSelection, noTimeCounterSubCycled, dictionary);



// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

//- Construct from components
noTimeCounterSubCycled::noTimeCounterSubCycled
(
    const polyMesh& mesh,
    dsmcCloud& cloud,
    const dictionary& dict
)
:
    collisionPartnerSelection(mesh, cloud, dict),
    propsDict_(dict.subDict(typeName + "Properties")),
    nSubCycles_(readLabel(propsDict_.lookup("nSubCycles"))),
    infoCounter_(0)
{}



// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

noTimeCounterSubCycled::~noTimeCounterSubCycled()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void noTimeCounterSubCycled::initialConfiguration()
{

}

void noTimeCounterSubCycled::collide()
{
    if (!cloud_.binaryCollision().active())
    {
        return;
    }

    // Temporary storage for subCells
    List<DynamicList<label> > subCells(8);

    scalar deltaT = cloud_.mesh().time().deltaTValue();

    label collisionCandidates = 0;

    label collisions = 0;

    const List<DynamicList<dsmcParcel*> > cellOccupancy = cloud_.cellOccupancy();

    const polyMesh& mesh = cloud_.mesh();
    
    for(label ii = 0; ii < nSubCycles_; ii++)
    {
        forAll(cellOccupancy, cellI)
        {
            const DynamicList<dsmcParcel*>& cellParcels(cellOccupancy[cellI]);
            const scalar& cellVolume = mesh.cellVolumes()[cellI];

            const label nC(cellParcels.size());

            if (nC > 1)
            {
                // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                // Assign particles to one of 8 Cartesian subCells

                // Clear temporary lists
                forAll(subCells, i)
                {
                    subCells[i].clear();
                }

                // Inverse addressing specifying which subCell a parcel is in
                List<label> whichSubCell(cellParcels.size());

                const point& cC = mesh.cellCentres()[cellI];

                forAll(cellParcels, i)
                {
                    const dsmcParcel& p = *cellParcels[i];

                    vector relPos = p.position() - cC;

                    label subCell =
                        pos(relPos.x()) + 2*pos(relPos.y()) + 4*pos(relPos.z());

                    subCells[subCell].append(i);

                    whichSubCell[i] = subCell;
                }

                // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

                scalar sigmaTcRMax = cloud_.sigmaTcRMax()[cellI];

                const scalar& RWF = cloud_.getRWF_cell(cellI);
                
                scalar selectedPairs =
                    cloud_.collisionSelectionRemainder()[cellI]
                    + 0.5*nC*(nC - 1)*cloud_.nParticle()*RWF*sigmaTcRMax*(deltaT/nSubCycles_)
                    /cellVolume;
                
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

                    const label nSC = subCellPs.size();

                    if (nSC > 1)
                    {
                        // If there are two or more particle in a subCell, choose
                        // another from the same cell.  If the same candidate is
                        // chosen, choose again. If two electrons are chosen,
                        // choose again.
                        
                        do
                        {
                            candidateQ = subCellPs[rndGen_.integer(0, nSC - 1)];

                        } while (candidateP == candidateQ);
                    }
                    else
                    {
                        // Select a possible second collision candidate from the
                        // whole cell.  If the same candidate is chosen, choose
                        // again. If two electrons are chosen, choose again.

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

                    label chargeP = -2;
                    label chargeQ = -2;

                    chargeP = cloud_.constProps(parcelP.typeId()).charge();
                    chargeQ = cloud_.constProps(parcelQ.typeId()).charge();
                    
                    //do not allow electron-electron collisions
                    
                    if(!(chargeP == -1 && chargeQ == -1))
                    {
                        scalar sigmaTcR = cloud_.binaryCollision().sigmaTcR
                        (
                            parcelP,
                            parcelQ
                        );
                        
//                         Pout << "sigmaTcR = " << sigmaTcR << endl;

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
        }
        
        cloud_.reBuildCellOccupancy();
    }

    reduce(collisions, sumOp<label>());

    reduce(collisionCandidates, sumOp<label>());

    cloud_.sigmaTcRMax().correctBoundaryConditions();
    
    infoCounter_++;
    
    if(infoCounter_ >= cloud_.nTerminalOutputs())
    {
        if (collisionCandidates)
        {
            Info<< "    Collisions                      = "
                << collisions << nl
//                 << "    Acceptance rate                 = "
//                 << scalar(collisions)/scalar(collisionCandidates) << nl
                << endl;
                
            infoCounter_ = 0;
        }
        else
        {
            Info<< "    No collisions" << endl;
            
            infoCounter_ = 0;
        }
    }
}

// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //



// * * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
