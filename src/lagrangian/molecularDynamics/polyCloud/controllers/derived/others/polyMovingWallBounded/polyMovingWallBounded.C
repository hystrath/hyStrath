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

#include "polyMovingWallBounded.H"
#include "addToRunTimeSelectionTable.H"
#include "IFstream.H"
#include "graph.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(polyMovingWallBounded, 0);
addToRunTimeSelectionTable(polyStateController, polyMovingWallBounded, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
polyMovingWallBounded::polyMovingWallBounded
(
    Time& t,
    polyMoleculeCloud& molCloud,
    const dictionary& dict
)
:
    polyStateController(t, /*mesh,*/ molCloud, dict),
    propsDict_(dict.subDict(typeName + "Properties")),
    wallMotionModel_(),
    molIds_()
{
	writeInTimeDir_ = false;
    writeInCase_ = false;

//     singleValueController() = true;

    wallMotionModel_ = autoPtr<wallMotion>
    (
       wallMotion::New(t, propsDict_)
    );

    molIds_.clear();

    selectIds ids
    (
        molCloud_.cP(),
        propsDict_
    );

    molIds_ = ids.molIds();

    vector startPoint = propsDict_.lookup("startPoint");
    vector endPoint = propsDict_.lookup("endPoint");

    bb_.resetBoundedBox(startPoint, endPoint);    
    
    readProperties();
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

polyMovingWallBounded::~polyMovingWallBounded()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void polyMovingWallBounded::initialConfiguration()
{
    {
        IDLList<polyMolecule>::iterator mol(molCloud_.begin());
        
        for (mol = molCloud_.begin(); mol != molCloud_.end(); ++mol)
        {
            if(findIndex(molIds_, mol().id()) != -1)
            {
                if(bb_.contains(mol().position()))
                {  
                    mol().special() = 0; // unfreeze
                }
            }
        }
    }    
    
}

void polyMovingWallBounded::controlBeforeVelocityI()
{
    wallMotionModel_->updateVelocity();
    velocity_ = wallMotionModel_->velocity();

//     Info << "polyMovingWallBounded: control velocity = " << velocity_ << endl;
        
    {
        IDLList<polyMolecule>::iterator mol(molCloud_.begin());
        
        for (mol = molCloud_.begin(); mol != molCloud_.end(); ++mol)
        {
            if(findIndex(molIds_, mol().id()) != -1)
            {
                if(bb_.contains(mol().position()))
                {  
                    mol().a() = vector::zero;
                    mol().v() = velocity_;                    
                }
            }
        }
    }   
}

void polyMovingWallBounded::controlBeforeMove()
{}

void polyMovingWallBounded::controlBeforeForces()
{}


void polyMovingWallBounded::controlDuringForces
(
    polyMolecule* molI,
    polyMolecule* molJ
)
{}

void polyMovingWallBounded::controlAfterForces()
{}

void polyMovingWallBounded::controlAfterVelocityII()
{}


void polyMovingWallBounded::calculateProperties()
{}

void polyMovingWallBounded::output
(
    const fileName& fixedPathName,
    const fileName& timePath
)
{}

void polyMovingWallBounded::updateProperties(const dictionary& newDict)
{
    //- the main controller properties should be updated first
    updateStateControllerProperties(newDict);

    propsDict_ = newDict.subDict(typeName + "Properties");

    readProperties();
}


void polyMovingWallBounded::readProperties()
{}

} // End namespace Foam

// ************************************************************************* //
