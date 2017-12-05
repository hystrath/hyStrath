/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2009-2010 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

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

\*---------------------------------------------------------------------------*/

#include "gravitationalController.H"
#include "addToRunTimeSelectionTable.H"
#include "fvc.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    
defineTypeNameAndDebug(gravitationalController, 0);

addToRunTimeSelectionTable(dsmcStateController, gravitationalController, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

gravitationalController::gravitationalController
(
    Time& t,
    dsmcCloud& cloud,
    const dictionary& dict
)
:
    dsmcStateController(t, cloud, dict),
    propsDict_(dict.subDict(typeName + "Properties")),
    model_(),
    nParcels_(0.0),
    acc_(vector::zero),
    nTimeSteps_(0.0),
    accCum_(vector::zero),
    nParcelsCum_(0.0)

{
    writeInTimeDir_ = false;
    writeInCase_ = false;
    singleValueController() = true;
    
    model_ = autoPtr<gravityForce>
    (
        gravityForce::New(t, propsDict_)
    );
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

gravitationalController::~gravitationalController()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void gravitationalController::initialConfiguration()
{
}

void gravitationalController::calculateProperties()
{
    vector acc = acc_;
    scalar nParcels = nParcels_;
    
    if(Pstream::parRun())
    {
        reduce(acc, sumOp<vector>());
        reduce(nParcels, sumOp<scalar>());        
    }

    accCum_ += acc;
    nParcelsCum_ += nParcels;
    nTimeSteps_ += 1.0;
    
    if(nParcels > 0)
    {    
        Info << nl << "gravitationalController " << nl
            << "Instantaneous: total acc = " << acc
            << ", total parcels = " << nParcels
            << ", acc per parcel = " << acc/nParcels
            << nl 
            << "Cumulative: acc per parcel = " <<  accCum_/nParcelsCum_
            << ", average acc per time-step = " << accCum_/nTimeSteps_
            << endl;
    }    
    
    acc_ = vector::zero;
    nParcels_ = 0.0;
}

void gravitationalController::controlParcelsBeforeMove()
{
    model_->updateForce();
        
    forAll(controlZone(), c)
    {
        const List<DynamicList<dsmcParcel*> >& cellOccupancy = cloud_.cellOccupancy();
        const label& cellI = controlZone()[c];
        const List<dsmcParcel*>& molsInCell = cellOccupancy[cellI];

        forAll(molsInCell, mIC)
        {
            dsmcParcel* p = molsInCell[mIC];

            vector acceleration = model_->force(p->position());
            
            p->U() += 0.5*acceleration*mesh_.time().deltaTValue();
            
            if(mag(acceleration) > SMALL)
            {
                acc_ += acceleration;
                nParcels_ += 1.0;
            }
        }
    }
}

void gravitationalController::output
(
    const fileName& fixedPathName,
    const fileName& timePath
)
{
    model_->write(fixedPathName, timePath);
}

void gravitationalController::controlParcelsBeforeCollisions()
{        
    forAll(controlZone(), c)
    {
        const List<DynamicList<dsmcParcel*> >& cellOccupancy = cloud_.cellOccupancy();
        const label& cellI = controlZone()[c];
        const List<dsmcParcel*>& molsInCell = cellOccupancy[cellI];
        forAll(molsInCell, mIC)
        {
            dsmcParcel* p = molsInCell[mIC];
            vector acceleration = model_->force(p->position());
            p->U() -= 0.5*acceleration*mesh_.time().deltaTValue();
        }
    }
}

void gravitationalController::controlParcelsAfterCollisions()
{
    forAll(controlZone(), c)
    {
        const List<DynamicList<dsmcParcel*> >& cellOccupancy = cloud_.cellOccupancy();
        const label& cellI = controlZone()[c];
        const List<dsmcParcel*>& molsInCell = cellOccupancy[cellI];
        forAll(molsInCell, mIC)
        {
            dsmcParcel* p = molsInCell[mIC];
            vector acceleration = model_->force(p->position());
            p->U() += acceleration*mesh_.time().deltaTValue();
        }
    }
}

void gravitationalController::updateProperties(const dictionary& newDict)
{
    updateStateControllerProperties(newDict);
    propsDict_ = newDict.subDict(typeName + "Properties");
    setProperties();
}

void gravitationalController::setProperties()
{

}


}

// ************************************************************************* //
