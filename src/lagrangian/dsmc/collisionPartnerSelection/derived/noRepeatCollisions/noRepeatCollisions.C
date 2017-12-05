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
    noRepeatCollisions

Description

\*----------------------------------------------------------------------------*/

#include "noRepeatCollisions.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(noRepeatCollisions, 0);

addToRunTimeSelectionTable(collisionPartnerSelection, noRepeatCollisions, dictionary);



// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

label noRepeatCollisions::pickFromCandidateList
(
    DynamicList<label>& candidatesInCell
)
{
//     Info << " list (before) " << candidatesInCell << endl;

    label entry = -1;
    label size = candidatesInCell.size();

    if(size > 0)
    {
        // choose a random number between 0 and the size of the candidateList size
        label randomIndex = rndGen_.integer(0, size - 1);
        entry = candidatesInCell[randomIndex];

//         Info<< "random index: " << randomIndex <<" entry " 
//             << entry << endl;

        // build a new list without the chosen entry

        DynamicList<label> newCandidates(0);
    
        forAll(candidatesInCell, i)
        {
            if(i != randomIndex)
            {
                newCandidates.append(candidatesInCell[i]);
            }
        }

        // transfer the new list
        candidatesInCell.transfer(newCandidates);
        candidatesInCell.shrink();

//         Info <<  " list (after) " << candidatesInCell << endl;
    }

    return entry;
}

void noRepeatCollisions::updateCandidateSubList
(
    const label& candidate,
    DynamicList<label>& candidatesInSubCell
)
{
//     Info << " updating sub list (before) " << candidatesInSubCell << endl;

    label newIndex = findIndex(candidatesInSubCell, candidate);

    DynamicList<label> newCandidates(0);

    forAll(candidatesInSubCell, i)
    {
        if(i != newIndex)
        {
            newCandidates.append(candidatesInSubCell[i]);
        }
    }

    // transfer the new list
    candidatesInSubCell.transfer(newCandidates);
    candidatesInSubCell.shrink();

//     Info <<  " list (after) " << candidatesInSubCell << endl;
}

label noRepeatCollisions::pickFromCandidateSubList
(
    DynamicList<label>& candidatesInCell,
    DynamicList<label>& candidatesInSubCell
)
{
//     Info << " list (before) " << candidatesInCell << endl;
//     Info << " sub list (before) " << candidatesInSubCell << endl;


    label entry = -1;
    label subCellSize = candidatesInSubCell.size();
    
    if(subCellSize > 0)
    {
        label randomIndex = rndGen_.integer(0, subCellSize - 1);
        entry = candidatesInSubCell[randomIndex];

//         Info<< "random index: " << randomIndex <<" entry " 
//             << entry << endl;

        DynamicList<label> newSubCellList(0);

        forAll(candidatesInSubCell, i)
        {
            if(i != randomIndex)
            {
                newSubCellList.append(candidatesInSubCell[i]);
            }
        }

        candidatesInSubCell.transfer(newSubCellList);
        candidatesInSubCell.shrink();

//         Info <<  " sub list (after) " << candidatesInSubCell << endl;

        label newIndex = findIndex(candidatesInCell, entry);

        DynamicList<label> newList(0);
    
        forAll(candidatesInCell, i)
        {
            if(i != newIndex)
            {
                newList.append(candidatesInCell[i]);
            }
        }

        candidatesInCell.transfer(newList);
        candidatesInCell.shrink();

//         Info <<  " list (after) " << candidatesInCell << endl;
    }

    return entry;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

//- Construct from components
noRepeatCollisions::noRepeatCollisions
(
    const polyMesh& mesh,
    dsmcCloud& cloud,
    const dictionary& dict
)
:
    collisionPartnerSelection(mesh, cloud, dict)
//     propsDict_(dict.subDict(typeName + "Properties"))
{}



// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

noRepeatCollisions::~noRepeatCollisions()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void noRepeatCollisions::initialConfiguration()
{

}

void noRepeatCollisions::collide()
{
    if (!cloud_.binaryCollision().active())
    {
        return;
    }

    // TEMP
//     label counter = 0;

    // Temporary storage for subCells
    List<DynamicList<label> > subCells(8);

    scalar deltaT = cloud_.mesh().time().deltaTValue();

    label collisionCandidates = 0;

    label collisions = 0;
	
	const List<DynamicList<dsmcParcel*> > cellOccupancy = cloud_.cellOccupancy();

    forAll(cellOccupancy, cellI)
    {
        const DynamicList<dsmcParcel*>& cellParcels(cellOccupancy[cellI]);

        label nC(cellParcels.size());

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

            const point& cC = mesh_.cellCentres()[cellI];

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

            scalar selectedPairs =
                cloud_.collisionSelectionRemainder()[cellI]
              + 0.5*nC*(nC - 1)*cloud_.nParticle()*sigmaTcRMax*deltaT
               /mesh_.cellVolumes()[cellI];

            label nCandidates(selectedPairs);

            cloud_.collisionSelectionRemainder()[cellI] = selectedPairs - nCandidates;

            collisionCandidates += nCandidates;

            //*********

            // list of candidates in cell
            DynamicList<label> candidateList(0);

            for (label c = 0; c < nC; c++)
            {
                candidateList.append(c);
            }

            candidateList.shrink();
//             counter++;

//             if(counter == 1)
//             {
//                 Info << nl << "cell: " << cellI << " start nCandidates: " << nCandidates << endl;

//                 Info << "candidateList: " << candidateList << endl;
//             }

            // list of candidates subcells
            List<DynamicList<label> > candidateSubList(subCells);

//             if(counter == 1)
//             {
//                 Info << "candidatesubcellList: " << candidateSubList << endl;
//             }
            //************

            for (label c = 0; c < nCandidates; c++)
            {

                // only-one candidate selection from subcells (OUR TRIAL VERSION)
                // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

                // Select the first collision candidate
                label candidateP = pickFromCandidateList(candidateList);

                // proceed with collision algorithm only if you find a candidate P.
                if(candidateP != -1)
                {
                    label sC = whichSubCell[candidateP];

                    DynamicList<label>& subCellPs = candidateSubList[sC];

                    updateCandidateSubList(candidateP, subCellPs);
                    
                    // Declare the second collision candidate
                    label candidateQ = pickFromCandidateSubList(candidateList, subCellPs);


                    // if you don't find candidateQ in the same subCell as candidateP
                    // try searching for it in the cell
                    if(candidateQ == -1)
                    {
                        candidateQ = pickFromCandidateList(candidateList);
    
                        if(candidateQ != -1)
                        {
                            sC = whichSubCell[candidateQ];

                            DynamicList<label>& subCellQs = candidateSubList[sC];

                            updateCandidateSubList(candidateQ, subCellQs);
                        }
                    }

                    // proceed with collision algorthm only if you find candidateQ 
                    if(candidateQ != -1)
                    {

//                         if(counter == 1)
//                         {
//                             Info << "candidateP: " << candidateP
//                                  << " candidateQ: " << candidateQ
//                                 << endl;
//                         }
                        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                        // uniform candidate selection procedure
        
        //                 // Select the first collision candidate
        //                 label candidateP = rndGen_.integer(0, nC-1);
        // 
        //                 // Select a possible second collision candidate
        //                 label candidateQ = rndGen_.integer(0, nC-1);
        // 
        //                 // If the same candidate is chosen, choose again
        //                 while (candidateP == candidateQ)
        //                 {
        //                     candidateQ = rndGen_.integer(0, nC-1);
        //                 }
        
        
        
        
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
                                        parcelQ,
                                        cellI
                                    );                                    
                                }
                                // if reaction unsuccessful use conventional collision model
                                if(cloud_.reactions().reactions()[rMId]->relax())
                                {
                                    cloud_.binaryCollision().collide
                                    (
                                        parcelP,
                                        parcelQ
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
                            
//                             Info << "Performed collision." << endl;
        
                            collisions++;
                        }
                    }
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
