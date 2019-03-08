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
    noTimeCounterTCC

Description

\*----------------------------------------------------------------------------*/

#include "noTimeCounterTCC.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(noTimeCounterTCC, 0);

addToRunTimeSelectionTable(collisionPartnerSelection, noTimeCounterTCC, dictionary);



// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

//- Construct from components
noTimeCounterTCC::noTimeCounterTCC
(
    const polyMesh& mesh,
    pdCloud& cloud,
    const dictionary& dict
)
:

    collisionPartnerSelection(mesh, cloud, dict),
    coeffDict_(dict.subDict(typeName + "Coeffs")),
    collisionCellOccupancy_(),
    collGroupList_(),
    ppC_(readScalar(coeffDict_.lookup("ParticlesPerCell"))),
    csL_(readScalar(coeffDict_.lookup("CellSizeLimit"))),
    //adaptInterval_(1),
    //intervalCounter_(0),
    debug_(false),
    firstRun_(true),
    pdColl_
    (
        IOobject
        (
            "pdColl",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero",dimless,0.0),
        zeroGradientFvPatchField<scalar>::typeName
    ),
    pdRndColl_
    (
        IOobject
        (
            "pdRndColl",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero",dimless,0.0),
        zeroGradientFvPatchField<scalar>::typeName
    ),
    pdRhoNColl_
    (
        IOobject
        (
            "pdRhoNColl",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero",dimless,0.0),
        zeroGradientFvPatchField<scalar>::typeName
    ),
    mfpReferenceTemperature_(273.0),
    typeIds_(),
    translationalT_
    (
        IOobject
        (
            "translationalT_NTC",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("0.0", dimTemperature, 0.0),
        zeroGradientFvPatchField<scalar>::typeName
    ),
    rhoN_
    (
        IOobject
        (
            "rhoN_NTC",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("0.0", dimless, 0.0),
        zeroGradientFvPatchField<scalar>::typeName
    ),
    rhoM_
    (
        IOobject
        (
            "rhoM_NTC",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("0.0", dimless/dimVolume, 0.0),
        zeroGradientFvPatchField<scalar>::typeName
    ),
    linearKE_
    (
        IOobject
        (
            "linearKE_NTC",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero",  dimensionSet(1, -1, -2, 0, 0), 0.0),
        zeroGradientFvPatchField<scalar>::typeName
    ),
    momentum_
    (
        IOobject
        (
            "momentum_NTC",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedVector
        (
            "zero",
            dimensionSet(1, -2, -1, 0, 0),
            vector::zero
        ),
        zeroGradientFvPatchField<vector>::typeName
    ),
    meanFreePath_
    (
        IOobject
        (
            "mfpVHS_ANTC",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero", dimLength, 0.0)
    ),
    nParcels_(),
    mfp_()
{

    //- initialise containers for VHS mean free path calculation
    typeIds_.setSize(cloud_.typeIdList().size(),-1);

    forAll(typeIds_,iD)
    {
        typeIds_[iD] = iD;
    }

    //- outer containers is type iD
    mfp_.setSize(typeIds_.size());
    nParcels_.setSize(typeIds_.size());

    //- set inner list size as mesh size
    forAll(mfp_, i)
    {
        nParcels_[i].setSize(mesh_.nCells());
        mfp_[i].setSize(mesh_.nCells());
    }

    if(coeffDict_.found("mfpReferenceTemperature"))
    {
        mfpReferenceTemperature_ = readScalar(coeffDict_.lookup("mfpReferenceTemperature"));
    }

    //if(coeffDict_.found("AdaptationInterval"))
    //{
    //    adaptInterval_ = readInt(coeffDict_.lookup("AdapationInterval"));
    //}

    if(coeffDict_.found("debug"))
    {
        debug_ = Switch(coeffDict_.lookup("debug"));
    }


}



// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

noTimeCounterTCC::~noTimeCounterTCC()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void noTimeCounterTCC::initialConfiguration()
{

}

void noTimeCounterTCC::resetFields()
{
    //- congolomerated mesh
    pdColl_        = dimensionedScalar("zero",dimless,0.0);
    pdRndColl_     = dimensionedScalar("zero",dimless,0.0);
    pdRhoNColl_    = dimensionedScalar("zero",dimless,0.0);

    forAll(collisionCellOccupancy_, cO)
    {
        collisionCellOccupancy_[cO].clear();
    }
    collGroupList_.clear();

    //- mean free path
    translationalT_ = dimensionedScalar("zero", dimTemperature,VSMALL);
    rhoN_           = dimensionedScalar("zero", dimless,VSMALL);
    rhoM_           = dimensionedScalar("zero", dimless/dimVolume,VSMALL);
    linearKE_       = dimensionedScalar("zero", dimensionSet(1, -1, -2, 0, 0), 0.0);
    momentum_       = dimensionedVector("zero", dimensionSet(1, -2, -1, 0, 0), vector::zero);

    meanFreePath_   = dimensionedScalar("zero",dimLength,0.0);

    //- reset inner list for mfp and mcr
    forAll(mfp_, iD)
    {
        mfp_[iD] = 0.0;
        nParcels_[iD] = 0.0;
    }

    //- reset other fields <-- can I just use meanFreePath_ = 0.0?
    forAll(meanFreePath_,cellI)
    {
        //- reset counters
        meanFreePath_[cellI] = 0.0;
    }
}


void noTimeCounterTCC::updateCollisionCellOccupancy()
{

    const fvMesh& mesh = cloud_.mesh();
    const List<DynamicList<pdParcel*> > cellOccupancy = cloud_.cellOccupancy();

    scalar familyN = 0.0;
    /** add in random parent selection **/
    //loop through mesh
    forAll(mesh_.cells(),cellI)
    {
        //start a new family?
        if(pdColl_[cellI] == 0.0)
        {
            //- reset particle counter
            scalar tP = 0;

            //- assign 0th parent a new family number
            scalar parent = 0;
            labelList familyList(1,cellI);

            familyN++;

            pdColl_[cellI] = familyN;

            //- cell dimension is measured relative to parent cell;
            scalar famDim = 0;
            scalar avemfp = meanFreePath_[cellI];

           //- add 0th parent particles to counter
            tP += cellOccupancy[cellI].size();

            while((tP<ppC_) && (famDim < avemfp/2) )
            {
                //- possible children for parent
                const labelList childList = mesh.cellCells()[familyList[parent]];
                //- maximum possible number of children
                int sizeC = childList.size();

                //- container list for children that are born
                labelList tempList(sizeC,0);
                //- counter for number of children born
                int births = 0;

                forAll(childList,cI)
                {
                    if(tP > ppC_)
                    {
                            //family is large enough :)
                            break;
                    }

                    //if cell is not already a families child
                    if(pdColl_[childList[cI]] == 0)
                    {
                        //- assign cell to parent's family
                        pdColl_[childList[cI]] = familyN;
                        //- add particles in cell to family counter
                        tP += cellOccupancy[childList[cI]].size();
                        //- record birth of child and which child it was
                        tempList[births] = childList[cI];
                        //- add to birth counter
                        births++;

                        //- check maximum cell dimension
                        scalar tempDim = mag(mesh_.C()[cellI]-mesh_.C()[childList[cI]])/2;
                        if(tempDim > famDim)
                        {
                            famDim = tempDim;
                        }
                    }
                }
                //check to see if there were any births
                if(births != 0)
                {
                    // strip non children from genList to get correct size
                    // zeros are at the end tempList so just iterate to size
                    // or births
                    labelList genList(births,0);
                    forAll(genList,cG)
                    {
                        genList[cG] = tempList[cG];
                    }

                    //append the new generation to the family list
                    familyList.append(genList);
                }

                //- get average mean free path from supercell
                scalar tempmfp = 0;
                forAll(familyList,member)
                {
                    if(meanFreePath_[familyList[member]] >= GREAT)
                    {
                        tempmfp += csL_;
                    }
                    else
                    {
                        tempmfp += meanFreePath_[familyList[member]];
                    }

                }
                avemfp = tempmfp/familyList.size();

                //move to the next parent in family list
                parent++;
                if(parent>=familyList.size())
                {
                    //couldn't fill out family, exit :(
                    break;
                }
            }

            forAll(familyList,member)
            {
                pdRhoNColl_[familyList[member]] = tP;
            }

            collGroupList_.append(familyList);
        }
    }

    collisionCellOccupancy_.resize(familyN);

    forAll(collGroupList_,family)
    {
        scalar rndNum = rndGen_.sample01<scalar>();
        labelList cellList = collGroupList_[family];
        scalar collVolume = 0.0;

        forAll(cellList,member)
        {
            //- randomise family groups for visualisation. wasn't done before because could not ensure unique ID
            const label& cellM = cellList[member];
            pdRndColl_[cellM] = rndNum;

            //- build collision cell occupancy list: outer list = super cell iD, inner list = particle iD
            forAll(cellOccupancy[cellM],partID)
            {
                collisionCellOccupancy_[family].append(cellOccupancy[cellM][partID]);
            }

            if(firstRun_)
            {
                collVolume += mesh_.cellVolumes()[cellM];
            }
        }

        if(firstRun_)
        {
            forAll(cellList,member)
            {
                const label& cellM = cellList[member];
                cloud_.collisionSelectionRemainder()[cellM] = rndNum*mesh_.cellVolumes()[cellM]/collVolume;
            }
        }
    }

    //- correct for processors and BCs
    pdRhoNColl_.correctBoundaryConditions();
    pdColl_.correctBoundaryConditions();

}

void noTimeCounterTCC::collide()
{
    if (!cloud_.binaryCollision().active())
    {
        return;
    }

    //- check if collision cells need to be updated
    //intervalCounter++;
    //if(intervalCounter == adaptInterval_)
    //{
         //- reset fields
        resetFields();
        measureMeanFreePath();
        updateCollisionCellOccupancy();
    //    intervalCounter = 0;
    //}

    //- collisionPartnerSelection metrics
    scalar deltaT = cloud_.mesh().time().deltaTValue();
    label collisionCandidates = 0;
    label collisions = 0;

    const List<DynamicList<pdParcel*> > cellOccupancy = cloud_.cellOccupancy();

    forAll(collisionCellOccupancy_, cO)
    {
        const DynamicList<pdParcel*>& cellParcels(collisionCellOccupancy_[cO]);
        labelList cellList = collGroupList_[cO];
        label nColl(cellParcels.size());                    //number of particles in collision cell

        if (nColl > 1)
        {
            /** Calculate number of collision pairs NTC **/

            //- get collision group properties
            scalar sigmaTcRMax = 0.0;
            scalar collVolume = 0.0;
            scalar collisionSelectionRemainder = 0.0;

            forAll(cellList,cG)
            {
                const label& cellG = cellList[cG];
                if(cloud_.sigmaTcRMax()[cellG]>sigmaTcRMax)
                {
                    sigmaTcRMax = cloud_.sigmaTcRMax()[cellG];
                }
                collisionSelectionRemainder += cloud_.collisionSelectionRemainder()[cellG];
                collVolume += mesh_.cellVolumes()[cellG];
            }

            //- calculate number of selected pairs
            scalar selectedPairs =
               collisionSelectionRemainder
               + 0.5*nColl*(nColl - 1)*cloud_.nParticle()*sigmaTcRMax*deltaT
               /collVolume;

            //- convert to integer value
            label nCandidates(selectedPairs);

            //- store collision selection remainder on mesh
            forAll(cellList,cG)
            {
                const label& cellG = cellList[cG];
                cloud_.collisionSelectionRemainder()[cellG] = (selectedPairs - nCandidates)*mesh_.cellVolumes()[cellG]/collVolume;
            }

            //- update total number of collision candidates
            collisionCandidates += nCandidates;

            //- attempt to collide particles
            for (label c = 0; c < nCandidates; c++)
            {
                // Supercell nearest partner collision selection
                // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                //- Select the first collision candidate
                //  label candidateP = rndGen_.position<label>(0, nColl - 1); OLD
                label candidateP = min(nColl-1,
                    label(rndGen_.sample01<scalar>()*nColl));

                //- Declare the second collision candidate
                label candidateQ = -1;
                label cellQ = -1;

                //- Declare candidateP cell iD
                const   label& cellP = cellParcels[candidateP]->cell();
                int     nSCP         = cellOccupancy[cellP].size();
                int     nSCQ         = 0;

                if (nSCP > 1)
                {
                    // If there are two or more particle in candidate's cell, choose
                    // another from the same cell.  If the same candidate is
                    // chosen, choose again.

                    do
                    {
                        //  candidateQ = rndGen_.position<label>(0, nSCP - 1); OLD
                        candidateQ = min(nSCP-1,
                            label(rndGen_.sample01<scalar>()*nSCP));

                        //transfer index of chosen particle in cellOccupancy to cellParcels list
                        forAll(cellParcels,cP)
                        {
                            if(cellParcels[cP] == cellOccupancy[cellP][candidateQ])
                            {
                                candidateQ = cP;
                                cellQ = cellP;
                                break;
                            }
                        }

                    //if chosen self, try again
                    } while (candidateP == candidateQ);
                }
                else
                {
                    // Select a possible second collision candidate from the
                    // nearest neighbour cell within the group with a candidate.
                    // If the same candidate is chosen, choose again.
                    labelList parentList(1,0);
                    parentList[0] = cellP;
                    label parent = 0;
                    label familyN = cO + 1;
                    while(nSCQ == 0)
                    {
                        //- possible children for parent
                        const labelList childList = cloud_.mesh().cellCells()[parentList[parent]];
                        forAll(childList,cN)
                        {
                            if(pdColl_[childList[cN]]==familyN)
                            {
                                if(cellOccupancy[childList[cN]].size() > 0)
                                {
                                    nSCQ = cellOccupancy[childList[cN]].size();
                                    cellQ = childList[cN];
                                    break;
                                }
                                parentList.append(childList[cN]);
                            }
                        }
                        parent++;

                        if(parent > parentList.size())
                        {
                            //couldn't find another particle
                            break;
                        }
                    }

                    if(cellQ != -1)
                    {
                        //  candidateQ = rndGen_.position<label>(0, nSCQ - 1); OLD
                        candidateQ = min(nSCQ-1,
                            label(rndGen_.sample01<scalar>()*nSCQ));// random number for cell

                        //transfer index of chosen particle in cellOccupancy to cellParcels list
                        forAll(cellParcels,cP)
                        {
                            if(cellParcels[cP] == cellOccupancy[cellQ][candidateQ])
                            {
                                candidateQ = cP;
                            }
                        }
                    }
                    else
                    {
                        do
                        {
                            //  candidateQ = rndGen_.position<label>(0, nColl - 1); OLD
                            candidateQ = min(nColl-1,
                                label(rndGen_.sample01<scalar>()*nColl));
                            cellQ = cellParcels[candidateQ]->cell();
                        } while (candidateP == candidateQ);

                        Info << "Warning: Could not find collision partner using TCC reverse stepping." << nl
                             << "Selecting random partner from within collision cell." << endl;
                    }
                }

                //de-referance candidates
                pdParcel& parcelP = *cellParcels[candidateP];
                pdParcel& parcelQ = *cellParcels[candidateQ];

                scalar sigmaTcR = cloud_.binaryCollision().sigmaTcR
                (
                    parcelP,
                    parcelQ
                );

                /** Update sigmaTcRMax in candidate subcells **/

                // Update the maximum value of sigmaTcR stored, but use the
                // initial value in the acceptance-rejection criteria because
                // the number of collision candidates selected was based on this
                if (sigmaTcR > sigmaTcRMax)
                {
                    forAll(cellList,cG)
                    {
                        const label& cellG = cellList[cG];
                        cloud_.sigmaTcRMax()[cellG] = sigmaTcR;
                    }
                }
                //Update collision group sigmaTcRMax
                //if (sigmaTcR > collGroupSigmaTcRMax_[cO])
                //{
                //    collGroupSigmaTcRMax_[cO] = sigmaTcR;
                //}
                /** Try and collide (acceptance - rejection) **/
                /****************************************************/

                if ((sigmaTcR/sigmaTcRMax) > rndGen_.sample01<scalar>())
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
                                parcelQ
                            );
                        }
                    }
                    else // if reaction model not found, use conventional collision model
                    {
                        cloud_.binaryCollision().collide
                        (
                            parcelP,
                            parcelQ
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
        Info<< "    Collision Cells                 = "
            << collGroupList_.size() << nl
            << "    Collisions                      = "
            << collisions << nl
            << "    Acceptance rate                 = "
            << scalar(collisions)/scalar(collisionCandidates) << nl
            << endl;
    }
    else
    {
        Info<< "    No collisions" << endl;
    }

    if(debug_)
    {
        rhoM_.write();
        rhoN_.write();
        linearKE_.write();
        momentum_.write();
        translationalT_.write();
    }
}

void noTimeCounterTCC::measureMeanFreePath()
{
    /** calculate VHS mean free path **/
    //- build list of number of type of parts in each cell subcell
    //- calculate required quantites
    forAllConstIter(pdCloud, cloud_, iter)
    {
        const pdParcel& p = iter();
        const label cellI = p.cell();
        label iD = findIndex(typeIds_, p.typeId());

        rhoN_[cellI]++;

        rhoM_[cellI] += cloud_.constProps(p.typeId()).mass();
        linearKE_[cellI] += 0.5*cloud_.constProps(p.typeId()).mass()*(p.U() & p.U());
        momentum_[cellI] += cloud_.constProps(p.typeId()).mass()*p.U();

        if(iD != -1)
        {
            nParcels_[iD][cellI] += 1.0;
        }
    }

    rhoN_.primitiveFieldRef() *= cloud_.nParticle()/mesh_.cellVolumes();
    rhoN_.correctBoundaryConditions();

    rhoM_.primitiveFieldRef() *= cloud_.nParticle()/mesh_.cellVolumes();
    rhoM_.correctBoundaryConditions();

    linearKE_.primitiveFieldRef() *= cloud_.nParticle()/mesh_.cellVolumes();
    linearKE_.correctBoundaryConditions();

    momentum_.primitiveFieldRef() *= cloud_.nParticle()/mesh_.cellVolumes();
    momentum_.correctBoundaryConditions();

    //- calculate mean free path for species and mixture
    forAll(translationalT_,cellI)
    {
        //- catch div0 error
        if(rhoN_[cellI] > SMALL)
        {
            vector cO = momentum_[cellI]/rhoM_[cellI];
            //- calculate translation temperature
            translationalT_[cellI] = 2.0/(3.0*physicoChemical::k.value()*rhoN_[cellI])
                                        *(linearKE_[cellI] - 0.5*rhoM_[cellI]*(cO & cO));

            //- calculate VHS mean free path for species
            forAll(mfp_, iD)
            {
                label qspec = 0;

                for (qspec=0; qspec < typeIds_.size(); qspec++)
                {
                    scalar dPQ = 0.5*(cloud_.constProps(typeIds_[iD]).d() + cloud_.constProps(typeIds_[qspec]).d());
                    scalar omegaPQ = 0.5*(cloud_.constProps(typeIds_[iD]).omega() + cloud_.constProps(typeIds_[qspec]).omega());
                    scalar massRatio = cloud_.constProps(typeIds_[iD]).mass()/cloud_.constProps(typeIds_[qspec]).mass();

                    if(nParcels_[qspec][cellI] > VSMALL && translationalT_[cellI] > VSMALL)
                    {
                        scalar nDensQ = (cloud_.nParticle()*nParcels_[qspec][cellI])/(mesh_.cellVolumes()[cellI]);

                        mfp_[iD][cellI] += (pi*dPQ*dPQ*nDensQ*pow(mfpReferenceTemperature_/translationalT_[cellI],omegaPQ-0.5)*sqrt(1.0+massRatio)); //TCC, eq (4.76)
                    }
                }

                if(mfp_[iD][cellI] > VSMALL)
                {
                    mfp_[iD][cellI] = 1.0/mfp_[iD][cellI];
                }
            }

            //- calculate mixture mean free path
            forAll(mfp_, iD)
            {
                if(rhoN_[cellI] > VSMALL)
                {
                    scalar nDensP = (cloud_.nParticle()*nParcels_[iD][cellI])/(mesh_.cellVolumes()[cellI]);

                    meanFreePath_[cellI] += mfp_[iD][cellI]*nDensP/rhoN_[cellI]; //TCC, eq (4.77)
                }
            }

            if(meanFreePath_[cellI] < VSMALL)
            {
                meanFreePath_[cellI] = GREAT;
            }
        }
        else
        {
            //- put user defined limit in
            if(meanFreePath_[cellI] < VSMALL)
            {
                meanFreePath_[cellI] = GREAT;
            }
        }
    }
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //



// * * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
