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

#include "forcingFunctionController.H"
#include "addToRunTimeSelectionTable.H"
#include "fvc.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(forcingFunctionController, 0);

    addToRunTimeSelectionTable(dsmcStateController, forcingFunctionController, dictionary);

    // * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

    forcingFunctionController::forcingFunctionController
    (
        Time& t,
        dsmcCloud& cloud,
        const dictionary& dict
    )
    :
    dsmcStateController(t, cloud, dict),
    propsDict_(dict.subDict(typeName + "Properties")),
    accelerationDirection_(propsDict_.lookup("accelerationDirection")),
    amplitude_(readScalar(propsDict_.lookup("amplitude"))),
    waveNumber_(readScalar(propsDict_.lookup("waveNumber"))),
    frequency_(readScalar(propsDict_.lookup("frequency"))),
    nTimeSteps_(0.0),
    currentTime_(0.0),
    initialXPositions_(),
    writeAcceleration_()
    {
        writeInTimeDir_ = false;
        writeInCase_ = false;
        singleValueController() = true;
    }

    // * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

    forcingFunctionController::~forcingFunctionController()
    {}

    // * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

    void forcingFunctionController::initialConfiguration()
    {
        initialXPositions_.setSize(cloud_.size());
        writeAcceleration_.setSize(mesh_.nCells());
    }

    void forcingFunctionController::calculateProperties()
    {
    }

    void forcingFunctionController::controlParcelsBeforeMove()
    {
        nTimeSteps_++;
        
        currentTime_ = mesh_.time().deltaTValue()*nTimeSteps_;
        
        forAll(controlZone(), c)
        {
            const List<DynamicList<dsmcParcel*> >& cellOccupancy = cloud_.cellOccupancy();
            const label& cellI = controlZone()[c];
            const List<dsmcParcel*>& molsInCell = cellOccupancy[cellI];

            forAll(molsInCell, mIC)
            {    
                dsmcParcel* p = molsInCell[mIC];
                
                label origID = p->origId();
                
                initialXPositions_[origID] = p->position().x();
                
                scalar accelerationMagnitude = amplitude_
                        *cos((waveNumber_*initialXPositions_[origID]) - (frequency_*currentTime_));
                        
                vector acceleration = accelerationDirection_*accelerationMagnitude;
                
                writeAcceleration_[c] = acceleration;
                
                p->U() += 0.5*acceleration*mesh_.time().deltaTValue();
            }
        }
    }

    void forcingFunctionController::output
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

    void forcingFunctionController::controlParcelsBeforeCollisions()
    {

        forAll(controlZone(), c)
        {
            const List<DynamicList<dsmcParcel*> >& cellOccupancy = cloud_.cellOccupancy();
            const label& cellI = controlZone()[c];
            const List<dsmcParcel*>& molsInCell = cellOccupancy[cellI];
            forAll(molsInCell, mIC)
            {
                dsmcParcel* p = molsInCell[mIC];
                
                label origID = p->origId();
                
                scalar accelerationMagnitude = amplitude_
                        *cos((waveNumber_*initialXPositions_[origID]) - (frequency_*currentTime_));
                        
                vector acceleration = accelerationDirection_*accelerationMagnitude;
                
                p->U() -= 0.5*acceleration*mesh_.time().deltaTValue();
            }
        }
    }

    void forcingFunctionController::controlParcelsAfterCollisions()
    {
        forAll(controlZone(), c)
        {
            const List<DynamicList<dsmcParcel*> >& cellOccupancy = cloud_.cellOccupancy();
            const label& cellI = controlZone()[c];
            const List<dsmcParcel*>& molsInCell = cellOccupancy[cellI];
            forAll(molsInCell, mIC)
            {
                dsmcParcel* p = molsInCell[mIC];
                
                label origID = p->origId();
                
                scalar accelerationMagnitude = amplitude_
                        *cos((waveNumber_*initialXPositions_[origID]) - (frequency_*currentTime_));
                        
                vector acceleration = accelerationDirection_*accelerationMagnitude;
                
                p->U() += acceleration*mesh_.time().deltaTValue();
            }
        }
    }

    void forcingFunctionController::updateProperties(const dictionary& newDict)
    {
        updateStateControllerProperties(newDict);
        propsDict_ = newDict.subDict(typeName + "Properties");
        setProperties();
    }

    void forcingFunctionController::setProperties()
    {      
        accelerationDirection_ = propsDict_.lookup("accelerationDirection");
        amplitude_ = readScalar(propsDict_.lookup("amplitude"));
        waveNumber_ = readScalar(propsDict_.lookup("waveNumber"));
        frequency_ = readScalar(propsDict_.lookup("frequency"));
    }
}

// ************************************************************************* //
