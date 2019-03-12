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

Control density in a user-defined zone.
This method provides good statistics, but loses resolution when inserting particles.

\*---------------------------------------------------------------------------*/

#include "densityZoneController.H"
#include "addToRunTimeSelectionTable.H"
#include "fvc.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(densityZoneController, 0);

addToRunTimeSelectionTable(dsmcStateController, densityZoneController, dictionary);



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
densityZoneController::densityZoneController
(
    Time& t,
    dsmcCloud& cloud,
    const dictionary& dict
)
:
    dsmcStateController(t, cloud, dict),
    propsDict_(dict.subDict(typeName + "Properties")),
    typeId_(-1),
    targetParcelDensity_(0.0),
    nMols_(controlZone().size(), 0),
    nMolsActualSum_(controlZone().size(), 0),
    avParcelDensity_(0.0),
    measuredParcels_(0.0),
    controlVolume_(0.0),
    volumeWeighting_()
{
    writeInTimeDir_ = false;
    writeInCase_ = false;

    singleValueController() = true;

}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

densityZoneController::~densityZoneController()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void densityZoneController::initialConfiguration()
{
    setProperties();

    forAll(controlZone(), c)
    {
        const label& cellI = controlZone()[c];
        controlVolume_ += mesh_.cellVolumes()[cellI];
    }


    if (Pstream::parRun())
    {
        volumeWeighting_.setSize(Pstream::nProcs(), 0.0);

        volumeWeighting_[Pstream::myProcNo()] = controlVolume_;

        for (int p = 0; p < Pstream::nProcs(); p++)
        {
            if(p != Pstream::myProcNo())
            {
                const int proc = p;
                {
                    OPstream toNeighbour(Pstream::commsTypes::blocking, proc);
                    toNeighbour << controlVolume_;
                }
            }
        }
    
        //- receiving
        for (int p = 0; p < Pstream::nProcs(); p++)
        {
            if(p != Pstream::myProcNo())
            {
                scalar controlVolumeProc;

                const int proc = p;
                {
                    IPstream fromNeighbour(Pstream::commsTypes::blocking, proc);
                    fromNeighbour >> controlVolumeProc;
                }

                volumeWeighting_[p] = controlVolumeProc;

                controlVolume_ += controlVolumeProc;
            }
        }

        forAll(volumeWeighting_, p)
        {
            volumeWeighting_[p] /= controlVolume_;
        }
    }
}

void densityZoneController::calculateProperties()
{

    if(time_.samplingTime())
    {
        const List<DynamicList<dsmcParcel*> >& cellOccupancy = cloud_.cellOccupancy();

        forAll(controlZone(), c)
        {
            const label& cellI = controlZone()[c];
            const List<dsmcParcel*>& parcelsInCell = cellOccupancy[cellI];

            forAll(parcelsInCell, pIC)
            {
                dsmcParcel* p = parcelsInCell[pIC];
                
                if(p->typeId() == typeId_)
                {
                    const scalar& RWF = cloud_.coordSystem().recalculateRWF(cellI);
                    
                    measuredParcels_ += RWF;
                }
            }
        }
    }

    if((time_.averagingTime()) && (controlZone().size() > 0))
    {
        if (Pstream::parRun())
        {
            reduce(measuredParcels_, sumOp<label>());
        }

        const scalar& nAvTimeSteps = time_.nAvTimeSteps().value();

        avParcelDensity_ = measuredParcels_/(controlVolume_*nAvTimeSteps);

        nMolsToControl();
        // reset
        measuredParcels_ = 0;
        nMolsActualSum_ = 0;
    }
}

void densityZoneController::nMolsToControl()
{
    nMols_ = scalar(0.0);

    scalar nMols = (targetParcelDensity_ - avParcelDensity_)*controlVolume_;

    // -total number of dsmcParcels to insert within entire control zone
    label nMolsRounded = 0;

    if(nMols > 0.0)
    {
        nMolsRounded = label(nMols + 0.5);
    }
    else if(nMols < 0.0)
    {
        nMolsRounded = label(nMols - 0.5);
    }

    //- Note (above): rounding off introduces an error

    // - number of dsmcParcels to insert within the decomposed control zone on each processor mesh
    label nMolsMesh = 0;

    if(Pstream::parRun())
    {
        //- distribute dsmcParcels based on the volume weighting of the control zone volume
        //  i.e. large zones insert more dsmcParcels than smaller zones

        if(nMolsRounded > 0.0)
        {
            nMolsMesh = label((scalar(nMolsRounded) * volumeWeighting_[Pstream::myProcNo()]) + 0.5);
        }
        else if(nMolsRounded < 0.0)
        {
            nMolsMesh = label((scalar(nMolsRounded) * volumeWeighting_[Pstream::myProcNo()]) - 0.5);
        }

        label nMolsMeshTotal = nMolsMesh;

        reduce(nMolsMeshTotal, sumOp<label>());

        // distribute residual dsmcParcels on the largest processor mesh

        scalar largerWeight = 0.0;

        label processorToInsert = -1;

        forAll(volumeWeighting_, p)
        {
            if(volumeWeighting_[p] > largerWeight)
            {
                largerWeight = volumeWeighting_[p];
                processorToInsert = p;
            }
        }

        if(processorToInsert == Pstream::myProcNo())
        {
            nMolsMesh += nMolsRounded - nMolsMeshTotal;
        }
    }
    else
    {
        nMolsMesh = nMolsRounded;
    }

    //- distribute dsmcParcels to cells on mesh

    if(controlZone().size() > 0)
    {
        label nMolsPerCell  = 0;

        if(nMolsMesh > 0.0)
        {
            nMolsPerCell = label((scalar(nMolsMesh)/scalar(controlZone().size()))/* + 0.5*/);
        }
        else if(nMolsMesh < 0.0)
        {
            nMolsPerCell = label((scalar(nMolsMesh)/scalar(controlZone().size()))/* - 0.5*/);
        }

        // More than one atomisticMolecule per cell
        if( mag(nMolsPerCell) >= 1 )
        {
            label nMolsCumul = 0;

            forAll(nMols_, c)
            {
                for (int n = 0; n < mag(nMolsPerCell); n++)
                {
                    if(mag(nMolsCumul) < mag(nMolsMesh))
                    {
                        if(nMolsMesh > 0)
                        {
                            nMols_[c] += 1;
                            nMolsCumul += 1;
                        }
                        else if(nMolsMesh < 0)
                        {
                            nMols_[c] -= 1;
                            nMolsCumul -= 1;
                        }
                    }
                }
            }

            //- residual dsmcParcels

            if( mag(nMolsCumul) < mag(nMolsMesh) )
            {
                label residualMols = nMolsMesh - nMolsCumul;
    
                DynamicList<label> cellsChosen(0);
            
                for (int n = 0; n < mag(residualMols); n++)
                {
                    label iter = 0;
                    bool foundCell = false;
        
                    while(!foundCell)
                    {
                        //label cellId = rndGen_.position<label>(0, controlZone().size()-1); OLD
                        label cellId = cloud_.randomLabel(0, controlZone().size()-1);
                        if( findIndex(cellsChosen, cellId) == -1)
                        {
                            cellsChosen.append(cellId);
                            foundCell = true;
                        }
                        else
                        {
                            iter++;
        
                            if(iter > controlZone().size())
                            {
                                cellsChosen.append(cellId);
                                foundCell = true;
                            }
                        }
                    }
                }
    
                cellsChosen.shrink();
    
                forAll(cellsChosen, c)
                {
                    if(mag(nMolsCumul) < mag(nMolsMesh))
                    {
                        if(nMolsMesh > 0)
                        {
                            nMols_[cellsChosen[c]] += 1;
                            nMolsCumul += 1;
                        }
                        else if(nMolsMesh < 0)
                        {
                            nMols_[cellsChosen[c]] -= 1;
                            nMolsCumul -= 1;
                        }
                    }
                }
            }
        }
        else //- residual dsmcParcels
        {
            DynamicList<label> cellsChosen(0);
        
            for (int n = 0; n < mag(nMolsMesh); n++)
            {
                label iter = 0;
                bool foundCell = false;
    
                while(!foundCell)
                {
                    //label cellId = rndGen_.position<label>(0, controlZone().size()-1); OLD
                    label cellId = cloud_.randomLabel(0, controlZone().size()-1);
    
                    if( findIndex(cellsChosen, cellId) == -1)
                    {
                        cellsChosen.append(cellId);
                        foundCell = true;
                    }
                    else
                    {
                        iter++;
    
                        if(iter > controlZone().size())
                        {
                            cellsChosen.append(cellId);
                            foundCell = true;
                        }
                    }
                }
            }

            cellsChosen.shrink();

            label nMolsCumul = 0;

            forAll(cellsChosen, c)
            {
                if(mag(nMolsCumul) < mag(nMolsMesh))
                {
                    if(nMolsMesh > 0)
                    {
                        nMols_[cellsChosen[c]] += 1;
                        nMolsCumul += 1;
                    }
                    else if(nMolsMesh < 0)
                    {
                        nMols_[cellsChosen[c]] -= 1;
                        nMolsCumul -= 1;
                    }
                }
            }
        }

        label totalMols = 0;

        forAll(nMols_, c)
        {
            totalMols += nMols_[c];
        }

        if(Pstream::parRun())
        {
            if(mag(totalMols) > 0)
            {
                Pout<< "sampled parcel density: " << avParcelDensity_ 
                    <<", nMols: " << nMols << ", mols per mesh: " <<  nMolsMesh 
                    << ", nMolsPerCell " << nMolsPerCell
                    << ", check on nMols : " << totalMols 
                    << endl;
            }
        }
        else
        {
            if(mag(totalMols) > 0)
            {
                Info<< "sampled parcel density: " << avParcelDensity_ 
                    <<", nMols: " << nMols << ", mols per mesh: " <<  nMolsMesh 
                    << ", nMolsPerCell " << nMolsPerCell
                    << ", check on nMols : " << totalMols 
                    << endl;
            }
        }
    }
}

void densityZoneController::controlParcelsBeforeMove()
{
    if((control_) && time_.controlTime())
    {
        Info << "Controlling density" << endl;

        const scalar& nControlSteps = time_.nControlSteps();

        labelField nMols(nMols_.size(), 0);
    
        forAll(nMols_, c)
        {
            if(nMols_[c] > 0)
            {
                nMols[c] = label((scalar(nMols_[c]) / nControlSteps) + 1.0);

                if((nMolsActualSum_[c] + nMols[c]) > nMols_[c])
                {
                    nMols[c] = nMols_[c] - nMolsActualSum_[c];
                }
            }
            else if (nMols_[c] < 0)
            {
                nMols[c] = label((scalar(nMols_[c]) / nControlSteps) - 1.0);
    
                if((nMolsActualSum_[c] + nMols[c]) < nMols_[c])
                {
                    nMols[c] = nMols_[c] - nMolsActualSum_[c];
                }
            }
        }

        forAll(nMols, c)
        {
            if(nMols[c] > 0) // insert parcels
            {
                insertParcels(nMols[c], c);
            }
            else if(nMols[c] < 0) // delete parcels
            {
                deleteParcels(nMols[c], c);
            }

            nMolsActualSum_[c] += nMols[c];
        }
    }
}

void densityZoneController::insertParcels(const label& nMols, const label& c)
{
    for (int n = 0; n < mag(nMols); n++)
    {
        const label& cellI = controlZone()[c];
    
        vector cC = mesh_.cellCentres()[cellI];
    
        // find the maximum distance between cell centre and cell vertices
        const labelList& cellPoints = mesh_.cellPoints()[cellI];
        scalar maxDistance = 0.0;
    
        forAll(cellPoints, cP)
        {
            const vector& vertexI = mesh_.points()[cellPoints[cP]];

            scalar vertexDist = mag(vertexI- cC);

            if(vertexDist > maxDistance)
            {
                maxDistance = vertexDist;
            }
        }
    
        // find a random point within the cell
        bool isPointInCell = false;

        vector p = vector::zero;

        while(!isPointInCell)
        {
            //- select a random direction
            vector randDirection = vector
            (
                rndGen_.GaussNormal<scalar>(),
                rndGen_.GaussNormal<scalar>(),
                rndGen_.GaussNormal<scalar>()
            );

            //- normalise the random vector (unit vector)
            randDirection /= mag(randDirection);
            
            p = randDirection*rndGen_.sample01<scalar>()*maxDistance + cC;

            if(mesh_.pointInCell(p, cellI))
            {
                isPointInCell = true;
            }
        }
    
        const dsmcParcel::constantProperties& cP = cloud_.constProps(typeId_);

        vector U = cloud_.equipartitionLinearVelocity
        (
            temperature_,
            cP.mass()
        );

        scalar ERot = cloud_.equipartitionRotationalEnergy
        (
            temperature_,
            cP.rotationalDegreesOfFreedom()
        );

        labelList vibLevel = cloud_.equipartitionVibrationalEnergyLevel
        (
            temperature_,
            cP.nVibrationalModes(),
            typeId_
        );
        
        label ELevel = cloud_.equipartitionElectronicLevel
        (
            temperature_,
            cP.electronicDegeneracyList(),
            cP.electronicEnergyList()
        );

        //thermal velocity + stream velocity = instantaneous velocity -
        // STREAM VELOCITY HERE IS USER DEFINED in the controllersDict.
        U += velocity_;

        label tetFace = -1;
        label tetPt = -1;

        mesh_.findTetFacePt
        (
            cellI,
            p,
            tetFace,
            tetPt
        );
        
        const scalar& RWF = cloud_.coordSystem().recalculateRWF(cellI);

        cloud_.addNewParcel
        (
            p,
            U,
            RWF,
            ERot,
            ELevel,
            cellI,
            tetFace,
            tetPt,
            typeId_,
            -1,
            0,
            vibLevel
        );
        
   // if parcel is inserted before the move function, the cell occupnacy need not be updated
    // only if the parcel is inserted after the buildCellOccupancy step, should the the 
    // cell occupancy be updated.
    }
} 

void densityZoneController::deleteParcels(const label& nMols, const label& c)
{
    for (int n = 0; n < mag(nMols); n++)
    {
        const List<DynamicList<dsmcParcel*> >& cellOccupancy = cloud_.cellOccupancy();
        const label& cellI = controlZone()[c];
        const List<dsmcParcel*>& molsInCell = cellOccupancy[cellI];
    
        if(molsInCell.size() > 0)
        {
            //label cellMolRemoveId = rndGen_.position<label>(0, molsInCell.size()-1);
            label cellMolRemoveId = cloud_.randomLabel(0, molsInCell.size()-1);
            dsmcParcel* delParcel = molsInCell[cellMolRemoveId];
            
            //- delete molecule from cellOccupancy (before deleting it from cloud)
            cloud_.removeParcelFromCellOccupancy(cellMolRemoveId, cellI);
            cloud_.deleteParticle(*delParcel);
        }
    }
}


void densityZoneController::output
(
    const fileName& fixedPathName,
    const fileName& timePath
)
{
    const Time& runTime = time_.time();

    if(runTime.outputTime())
    {

    }
}
void densityZoneController::controlParcelsBeforeCollisions()
{}

void densityZoneController::controlParcelsAfterCollisions()
{}

void densityZoneController::updateProperties(const dictionary& newDict)
{
    //- the main controller properties should be updated first
    updateStateControllerProperties(newDict);

    propsDict_ = newDict.subDict(typeName + "Properties");

    setProperties();
}

void densityZoneController::setProperties()
{
    word typeIdName = propsDict_.lookup("typeId");
    typeId_ = findIndex(cloud_.typeIdList(), typeIdName);

    if(typeId_ == -1)
    {
        FatalErrorIn("densityZoneController::densityZoneController()")
            << "Cannot find typeid: " << typeIdName << nl << "in: "
            << time_.time().system()/"controllersDict"
            << exit(FatalError);
    }

    velocity_ = propsDict_.lookup("velocity");
    temperature_ = readScalar(propsDict_.lookup("temperature"));

    targetParcelDensity_ = readScalar(propsDict_.lookup("massDensity"));

    scalar mass = cloud_.constProps(typeId_).mass();

    targetParcelDensity_ /= (mass*cloud_.nParticle());
}


}
// End namespace Foam

// ************************************************************************* //
