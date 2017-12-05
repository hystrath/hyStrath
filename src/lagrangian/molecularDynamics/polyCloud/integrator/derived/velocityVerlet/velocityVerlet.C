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

#include "velocityVerlet.H"
#include "addToRunTimeSelectionTable.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(velocityVerlet, 0);

addToRunTimeSelectionTable(polyIntegrator, velocityVerlet, dictionary);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
velocityVerlet::velocityVerlet
(
    Time& t,
    polyMoleculeCloud& molCloud,
    const dictionary& dict
)
:
    polyIntegrator(t, molCloud, dict)
/*    propsDict_(dict.subDict(typeName + "Properties")),*/    
  
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

velocityVerlet::~velocityVerlet()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void velocityVerlet::init()
{
    
}

void velocityVerlet::evolve()
{
    molCloud_.controlBeforeVelocity();
    updateVelocity(mesh_.time().deltaT().value());
    molCloud_.controlBeforeMove();
    molCloud_.move();
    molCloud_.controlAfterMove();
    molCloud_.buildCellOccupancy();
    molCloud_.controlBeforeForces();
    molCloud_.clearLagrangianFields();
    molCloud_.calculateForce();
    molCloud_.updateAcceleration();
    molCloud_.controlAfterForces();
    updateVelocity(mesh_.time().deltaT().value());
    molCloud_.controlAfterVelocity();
    molCloud_.postTimeStep();
}

void velocityVerlet::updateVelocity(const scalar& trackTime)
{
    IDLList<polyMolecule>::iterator mol(molCloud_.begin());

    for (mol = molCloud_.begin(); mol != molCloud_.end(); ++mol)
    {
        if(!mol().frozen())
        {
            mol().updateHalfVelocity(molCloud_.cP(), trackTime);
        }
    }
}


} // End namespace Foam

// ************************************************************************* //
