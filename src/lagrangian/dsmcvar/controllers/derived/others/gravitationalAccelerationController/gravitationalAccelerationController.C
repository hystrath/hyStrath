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

#include "gravitationalAccelerationController.H"
#include "addToRunTimeSelectionTable.H"
#include "fvc.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(gravitationalAccelerationController, 0);

    addToRunTimeSelectionTable(dsmcStateController, gravitationalAccelerationController, dictionary);

    // * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

    gravitationalAccelerationController::gravitationalAccelerationController
    (
        Time& t,
        dsmcCloud& cloud,
        const dictionary& dict
    )
    :
    dsmcStateController(t, cloud, dict),
    propsDict_(dict.subDict(typeName + "Properties")),
    acceleration_(propsDict_.lookup("acceleration"))
    {
        writeInTimeDir_ = false;
        writeInCase_ = false;
        singleValueController() = true;
    }

    // * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

    gravitationalAccelerationController::~gravitationalAccelerationController()
    {}

    // * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

    void gravitationalAccelerationController::initialConfiguration()
    {
    }

    void gravitationalAccelerationController::calculateProperties()
    {
    }

    void gravitationalAccelerationController::controlParcelsBeforeMove()
    {
        forAll(controlZone(), c)
        {
            const List<DynamicList<dsmcParcel*> >& cellOccupancy = cloud_.cellOccupancy();
            const label& cellI = controlZone()[c];
            const List<dsmcParcel*>& molsInCell = cellOccupancy[cellI];

            forAll(molsInCell, mIC)
            {
                dsmcParcel* p = molsInCell[mIC];
                p->U() += 0.5*acceleration_*mesh_.time().deltaTValue();
            }
        }
    }

    void gravitationalAccelerationController::output
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

    void gravitationalAccelerationController::controlParcelsBeforeCollisions()
    {
        forAll(controlZone(), c)
        {
            const List<DynamicList<dsmcParcel*> >& cellOccupancy = cloud_.cellOccupancy();
            const label& cellI = controlZone()[c];
            const List<dsmcParcel*>& molsInCell = cellOccupancy[cellI];
            forAll(molsInCell, mIC)
            {
                dsmcParcel* p = molsInCell[mIC];
                p->U() -= 0.5*acceleration_*mesh_.time().deltaTValue();
            }
        }
    }

    void gravitationalAccelerationController::controlParcelsAfterCollisions()
    {
        forAll(controlZone(), c)
        {
            const List<DynamicList<dsmcParcel*> >& cellOccupancy = cloud_.cellOccupancy();
            const label& cellI = controlZone()[c];
            const List<dsmcParcel*>& molsInCell = cellOccupancy[cellI];
            forAll(molsInCell, mIC)
            {
                dsmcParcel* p = molsInCell[mIC];
                p->U() += acceleration_*mesh_.time().deltaTValue();
            }
        }
    }

    void gravitationalAccelerationController::updateProperties(const dictionary& newDict)
    {
        updateStateControllerProperties(newDict);
        propsDict_ = newDict.subDict(typeName + "Properties");
        setProperties();
    }

    void gravitationalAccelerationController::setProperties()
    {
        acceleration_ = propsDict_.lookup("acceleration");
    }
}

// ************************************************************************* //
