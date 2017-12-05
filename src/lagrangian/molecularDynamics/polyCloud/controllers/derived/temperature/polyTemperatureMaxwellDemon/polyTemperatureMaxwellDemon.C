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

#include "polyTemperatureMaxwellDemon.H"
#include "addToRunTimeSelectionTable.H"
#include "IFstream.H"
#include "graph.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(polyTemperatureMaxwellDemon, 0);

addToRunTimeSelectionTable(polyStateController, polyTemperatureMaxwellDemon, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
polyTemperatureMaxwellDemon::polyTemperatureMaxwellDemon
(
    Time& t,
    polyMoleculeCloud& molCloud,
    const dictionary& dict
)
:
    polyStateController(t, molCloud, dict),
    propsDict_(dict.subDict(typeName + "Properties")),
    tauT_(readScalar(propsDict_.lookup("tauT"))),
    p_(0.0),
    molIds_()
{
    writeInTimeDir_ = false;
    writeInCase_ = false;

    velocities_.setSize(controlZone().size(), vector::zero);

    velocity_ = propsDict_.lookup("velocity");

    forAll(velocities_, c)
    {
        velocities_[c] = velocity_;
    }

//     fieldController() = true;

    temperature_ = readScalar(propsDict_.lookup("temperature"));

    temperatures_.setSize(controlZone().size(), 0.0);

    forAll(temperatures_, c)
    {
        temperatures_[c] = temperature_;
    }

    p_ = 1.0 - exp(-(time_.deltaT().value()/tauT_));

    Info << "probability of collisions:  " << p_ << endl;

    molIds_.clear();

    selectIds ids
    (
        molCloud_.cP(),
        propsDict_
    );

    molIds_ = ids.molIds();
    
    Info << "mol ids = " << molIds_ << endl;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

polyTemperatureMaxwellDemon::~polyTemperatureMaxwellDemon()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void polyTemperatureMaxwellDemon::initialConfiguration()
{}

void polyTemperatureMaxwellDemon::controlBeforeVelocityI()
{}

void polyTemperatureMaxwellDemon::controlBeforeMove()
{}

void polyTemperatureMaxwellDemon::controlBeforeForces()
{}

void polyTemperatureMaxwellDemon::controlDuringForces
(
    polyMolecule* molI,
    polyMolecule* molJ
)
{}

void polyTemperatureMaxwellDemon::controlAfterForces()
{}

void polyTemperatureMaxwellDemon::controlAfterVelocityII()
{
    if(control_)
    {
        Info << "polyTemperatureMaxwellDemon: control" << endl;

        forAll(controlZone(), c)
        {
            const label& cellI = controlZone()[c];
    
            const List<polyMolecule*>& molsInCell = molCloud_.cellOccupancy()[cellI];

            forAll(molsInCell, m)
            {
                polyMolecule* molI = molsInCell[m];

                if(findIndex(molIds_, molI->id()) != -1)
                {
                    if(molCloud_.rndGen().sample01<scalar>() <= p_)
                    {
                        const scalar& massI = molCloud_.cP().mass(molI->id());

                        scalar sigma = sqrt (temperatures_[c]*molCloud_.redUnits().kB() / massI);
    
                        vector molVel
                        (
                            sigma*molCloud_.rndGen().GaussNormalMD<scalar>(),
                            sigma*molCloud_.rndGen().GaussNormalMD<scalar>(),
                            sigma*molCloud_.rndGen().GaussNormalMD<scalar>()
                        );
    
                        molI->v() = velocities_[c] + molVel;
                    }
                }
            }
        }
    }
}


void polyTemperatureMaxwellDemon::calculateProperties()
{}

void polyTemperatureMaxwellDemon::output
(
    const fileName& fixedPathName, 
    const fileName& timePath
)
{}

void polyTemperatureMaxwellDemon::updateProperties(const dictionary& newDict)
{
    //- the main controller properties should be updated first
    updateStateControllerProperties(newDict);

    propsDict_ = newDict.subDict(typeName + "Properties");
/*
    if (propsDict_.found("tauT"))
    {
        tauT_ = readScalar(propsDict_.lookup("tauT"));
    }

    p_ = 1.0 - exp(-(time_.deltaT().value()/tauT_));

    Info << "probability of collisions:  " << p_ << endl;

    if (readStateFromFile_)
    {
        velocity_ = propsDict_.lookup("velocity");

        forAll(velocities_, c)
        {
            velocities_[c] = velocity_;
        }

        temperature_ = readScalar(propsDict_.lookup("temperature"));

        forAll(temperatures_, c)
        {
            temperatures_[c] = temperature_;
        }
    }*/
}

} // End namespace Foam

// ************************************************************************* //
