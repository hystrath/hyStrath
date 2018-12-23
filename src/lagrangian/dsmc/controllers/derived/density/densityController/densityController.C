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

#include "densityController.H"
#include "addToRunTimeSelectionTable.H"
#include "fvc.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(densityController, 0);

addToRunTimeSelectionTable(dsmcStateController, densityController, dictionary);



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
densityController::densityController
(
    Time& t,
    dsmcCloud& cloud,
    const dictionary& dict
)
:
    dsmcStateController(t, cloud, dict),
    propsDict_(dict.subDict(typeName + "Properties")),
    typeId_(-1),
    density_(readScalar(propsDict_.lookup("massDensity"))),
    velocity_(propsDict_.lookup("velocity")),
    temperature_(readScalar(propsDict_.lookup("temperature"))),
    densities_(controlZone().size(), 0.0),
    velocities_(controlZone().size(), vector::zero),
    temperatures_(controlZone().size(), 0.0),
    nMolsToControl_(controlZone().size(), 0.0),
    parcels_(controlZone().size(), 0.0),
    measuredDensity_(controlZone().size(), 0.0),
    residual_(controlZone().size(), 0.0),
    residualSum_(controlZone().size(), 0.0)

{

}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

densityController::~densityController()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void densityController::initialConfiguration()
{
	setProperties();
}

void densityController::calculateProperties()
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
                    
                    parcels_[c] += RWF;
                }
            }
        }
    }

    if((time_.averagingTime()) && (controlZone().size() > 0))
    {
        Info << "Averaging" << endl;
        const scalar& nAvTimeSteps = time_.nAvTimeSteps().value();
        Info << "nAvTimeSteps = " << nAvTimeSteps << endl;

        forAll(measuredDensity_, c)
        {
            measuredDensity_[c] = parcels_[c]/nAvTimeSteps;
            Info << "measuredDensity_[c] = " << measuredDensity_[c] << endl;
            Info << "densities_[c] = " << densities_[c] << endl;
            scalar deltaNParcels = densities_[c] - measuredDensity_[c];
            Info << "deltaNParcels = " << deltaNParcels << endl;
            scalar integerNParcels = 0.0;

            if(deltaNParcels > 0.0)
            {
                integerNParcels = label(deltaNParcels + 0.5);
            }
            else if (deltaNParcels < 0.0)
            {
                integerNParcels = label(deltaNParcels - 0.5);
            }

            nMolsToControl_[c] = integerNParcels;
            residual_[c] = deltaNParcels - nMolsToControl_[c];
            parcels_[c] = 0.0;
            Info << "Number of parcels to control = " << integerNParcels << endl;
            Info << "Number of residual parcels = " << residual_[c] << endl;
        }
    }
}


void densityController::controlParcelsBeforeMove()
{
    if((control_) && time_.controlTime())
    {
        Info << "Controlling density" << endl;

        labelField nMolsToControl(nMolsToControl_.size(), 0);
        scalarField residualSum(residualSum_.size(), 0);

        forAll(nMolsToControl_, c)
        {
            if(nMolsToControl_[c] > 0) // insert parcels
            {
                insertParcelWithinDSMC(c);
            }
            else if(nMolsToControl_[c] < 0) // delete parcels
            {
                deleteParcelFromDSMC(c);
            }
        }

        forAll(residualSum_, c)
        {
            residualSum_[c] += residual_[c];

            if(residualSum_[c] > 1) // insert residual parcel
            {
                insertParcelWithinDSMC(c);
                residualSum_[c] =- 1;
            }
            else if(residualSum_[c] < -1) // delete residual parcel
            {
                deleteParcelFromDSMC(c);
                residualSum_[c] =+ 1;
            }
        }
        
        nMolsToControl_ = 0;
    }
}



void densityController::deleteParcelFromDSMC(const label& c)
{
    const List<DynamicList<dsmcParcel*> >& cellOccupancy = cloud_.cellOccupancy();
    const label& cellI = controlZone()[c];
    const List<dsmcParcel*>& molsInCell = cellOccupancy[cellI];

    if(molsInCell.size() > 0)
    {
        //label cellMolRemoveId = rndGen_.position<label>(0, molsInCell.size()-1); OLD
        label cellMolRemoveId = rndGen_.sample01<scalar>()*molsInCell.size();
        dsmcParcel* delParcel = molsInCell[cellMolRemoveId];

        //- delete molecule from cellOccupancy (before deleting it from cloud)
        cloud_.removeParcelFromCellOccupancy(cellMolRemoveId, cellI);
        cloud_.deleteParticle(*delParcel);
    }

}

void densityController::insertParcelWithinDSMC(const label& c)
{
    const label& cellI = controlZone()[c];
    const scalar& cellTemperature = temperatures_[c];

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
        cP.electronicEnergyList(),
        typeId_
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
} 

void densityController::controlParcelsBeforeCollisions()
{
}

void densityController::controlParcelsAfterCollisions()
{
}

void densityController::output
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

void densityController::updateProperties(const dictionary& newDict)
{
    //- the main controller properties should be updated first
    updateStateControllerProperties(newDict);

    propsDict_ = newDict.subDict(typeName + "Properties");

    if(readStateFromFile_)
    {
        density_ = readScalar(propsDict_.lookup("massDensity"));
        velocity_ = propsDict_.lookup("velocity");
        temperature_ = readScalar(propsDict_.lookup("temperature"));

        density_ /= (cloud_.constProps(typeId_).mass() * cloud_.nParticle());

        const scalarField& volField = mesh_.cellVolumes();
        
        forAll(densities_, c)
        {
            const label& cellI = controlZone()[c];
            densities_[c] = density_*volField[cellI];
        }

        forAll(velocities_, c)
        {
            velocities_[c] = velocity_;
        }

        forAll(temperatures_, c)
        {
            temperatures_[c] = temperature_;
        }
    }
}

void densityController::setProperties()
{
    // search for mol id
    const List<word>& idList(cloud_.typeIdList());
    const word typeId = propsDict_.lookup("typeId");
    typeId_ = findIndex(cloud_.typeIdList(), typeId);

    if(typeId_ == -1)
    {
        FatalErrorIn("densityController::densityController()")
        << "Cannot find typeId: " << typeId << nl << "in: "
        << time_.time().system()/"controllersDict"
        << exit(FatalError);
    }

    Info << "typeId_ = " << typeId_ << endl;
    density_ /= (cloud_.constProps(typeId_).mass() * cloud_.nParticle()); // = number density
    const scalarField& volField = mesh_.cellVolumes();

    forAll(densities_, c)
    {
        const label& cellI = controlZone()[c];
        densities_[c] = density_*volField[cellI];
        Info << "densities_[c] = " << densities_[c] << endl;
    }

    forAll(velocities_, c)
    {
        velocities_[c] = velocity_;
    }

    forAll(temperatures_, c)
    {
        temperatures_[c] = temperature_;
    }
}

}
// End namespace Foam

// ************************************************************************* //
