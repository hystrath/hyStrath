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

#include "polyForce.H"
#include "addToRunTimeSelectionTable.H"
#include "IFstream.H"
#include "graph.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(polyForce, 0);
addToRunTimeSelectionTable(polyStateController, polyForce, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
polyForce::polyForce
(
    Time& t,
    polyMoleculeCloud& molCloud,
    const dictionary& dict
)
:
    polyStateController(t, molCloud, dict),
    propsDict_(dict.subDict(typeName + "Properties")),
    model_(),
    molIds_(),
    nTimeSteps_(0.0),
    force_(vector::zero)
{

    writeInTimeDir_ = true;
    writeInCase_ = true;

//     singleValueController() = true;

    model_ = autoPtr<gravityForce>
    (
        gravityForce::New(t, propsDict_)
    );

    molIds_.clear();

    selectIds ids
    (
        molCloud_.cP(),
        propsDict_
    );

    molIds_ = ids.molIds();

    readProperties();
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

polyForce::~polyForce()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void polyForce::initialConfiguration()
{}

void polyForce::controlBeforeVelocityI()
{}

void polyForce::controlBeforeMove()
{}

void polyForce::controlBeforeForces()
{}

void polyForce::controlDuringForces
(
    polyMolecule* molI,
    polyMolecule* molJ
)
{}

void polyForce::controlAfterForces()
{
    //set target velocity of wall
    model_->updateForce();

    // - if control switch is on
    if(control_) 
    {
        Info << "polyForce: control" << endl;

        forAll(controlZone(), c)
        {
            const label& cellI = controlZone()[c];
    
            const List<polyMolecule*>& molsInCell = molCloud_.cellOccupancy()[cellI];
    
            forAll(molsInCell, m)
            {
                polyMolecule* molI = molsInCell[m];
    
                if(findIndex(molIds_, molI->id()) != -1)
                {
                    const scalar& massI = molCloud_.cP().mass(molI->id());

                    vector force = model_->force(molI->position());

                    molI->a() += force/massI;

                    force_ += force;
                }
            }
        }

        nTimeSteps_ += 1.0;
    }
}

void polyForce::controlAfterVelocityII()
{}

void polyForce::calculateProperties()
{}

void polyForce::output
(
    const fileName& fixedPathName,
    const fileName& timePath
)
{
    model_->write(fixedPathName, timePath);
}

void polyForce::updateProperties(const dictionary& newDict)
{
    //- the main controller properties should be updated first
    updateStateControllerProperties(newDict);

    propsDict_ = newDict.subDict(typeName + "Properties");

    model_->updateProperties(propsDict_);

    readProperties();
}

void polyForce::readProperties()
{}

} // End namespace Foam

// ************************************************************************* //
