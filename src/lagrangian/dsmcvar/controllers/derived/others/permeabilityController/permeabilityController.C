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

#include "permeabilityController.H"
#include "addToRunTimeSelectionTable.H"
#include "fvc.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(permeabilityController, 0);

    addToRunTimeSelectionTable(dsmcStateController, permeabilityController, dictionary);

    // * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

    permeabilityController::permeabilityController
    (
        Time& t,
        dsmcCloud& cloud,
        const dictionary& dict
    )
    :
    dsmcStateController(t, cloud, dict),
    propsDict_(dict.subDict(typeName + "Properties")),
    porosity_(readScalar(propsDict_.lookup("porosity"))),
    temperature_(readScalar(propsDict_.lookup("surfaceTemperature")))
    {
        writeInTimeDir_ = false;
        writeInCase_ = false;
        singleValueController() = true;
    }

    // * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

    permeabilityController::~permeabilityController()
    {}

    // * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

    void permeabilityController::initialConfiguration()
    {
    }

    void permeabilityController::calculateProperties()
    {
    }

    void permeabilityController::controlParcelsBeforeMove()
    {

    }

    void permeabilityController::output
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

    void permeabilityController::controlParcelsBeforeCollisions()
    {
    }

    void permeabilityController::controlParcelsAfterCollisions()
    {
        forAll(controlZone(), c)
        {
            const List<DynamicList<dsmcParcel*> >& cellOccupancy = cloud_.cellOccupancy();
            const label& cellI = controlZone()[c];
            const List<dsmcParcel*>& molsInCell = cellOccupancy[cellI];
            Random& rndGen(cloud_.rndGen());
//             const scalar timeStepRatio = mesh_.time().deltaTValue()/1.0e-12;
            
            forAll(molsInCell, mIC)
            {
                dsmcParcel* p = molsInCell[mIC];
                
                if(rndGen.scalar01() > porosity_)
                {
                    scalar mass = cloud_.constProps(p->typeId()).mass();
                    
//                     Info << "Velocity before control = " << p->U() << endl;
                    
                    p->U() = sqrt(physicoChemical::k.value()*temperature_/mass)
                            *vector
                            (
                                rndGen_.GaussNormal(),
                                rndGen_.GaussNormal(),
                                rndGen_.GaussNormal()
                            );
                            
//                     Info << "Velocity after control = " << p->U() << endl;
                }
            }
        }
    }

    void permeabilityController::updateProperties(const dictionary& newDict)
    {
        updateStateControllerProperties(newDict);
        propsDict_ = newDict.subDict(typeName + "Properties");
        setProperties();
    }

    void permeabilityController::setProperties()
    {
        porosity_ = readScalar(propsDict_.lookup("porosity"));
        temperature_ = readScalar(propsDict_.lookup("surfaceTemperature"));
    }
}

// ************************************************************************* //
