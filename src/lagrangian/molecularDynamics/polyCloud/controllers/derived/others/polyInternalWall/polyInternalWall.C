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

#include "polyInternalWall.H"
#include "addToRunTimeSelectionTable.H"
#include "IFstream.H"
#include "graph.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(polyInternalWall, 0);

addToRunTimeSelectionTable(polyStateController, polyInternalWall, dictionary);






// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
polyInternalWall::polyInternalWall
(
    Time& t,
//     const polyMesh& mesh,
    polyMoleculeCloud& molCloud,
    const dictionary& dict
)
:
    polyStateController(t, /*mesh,*/ molCloud, dict),
    propsDict_(dict.subDict(typeName + "Properties")),
    
    molIds_()    
{

    writeInTimeDir_ = true;
    writeInCase_ = true;

//     singleValueController() = true;

    molIds_.clear();

    selectIds ids
    (
        molCloud_.cP(),
        propsDict_
    );

    molIds_ = ids.molIds();

    readProperties();


    startPoint_ = propsDict_.lookup("planeMidPoint");
    normal_ = propsDict_.lookup("planeNormal");
    
    normal_ /= mag(normal_);
    
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

polyInternalWall::~polyInternalWall()
{}



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void polyInternalWall::initialConfiguration()
{}

void polyInternalWall::controlBeforeVelocityI()
{}

void polyInternalWall::controlBeforeMove()
{}

void polyInternalWall::controlBeforeForces()
{}


void polyInternalWall::controlDuringForces
(
    polyMolecule* molI,
    polyMolecule* molJ
)
{}

void polyInternalWall::controlAfterForces()
{
    // - if control switch is on
    if(control_)
    {
        Info << "polyInternalWall: control" << endl;

        IDLList<polyMolecule>::iterator mol(molCloud_.begin());

        for (mol = molCloud_.begin(); mol != molCloud_.end(); ++mol)
        {
            if(findIndex(molIds_, mol().id()) != -1)
            {
                vector rSI = mol().position() - startPoint_;
                
                if((rSI & normal_) > 0)
                {
                    scalar vn = mol().v() & normal_;
                    
                    if(vn > 0)
                    {
                        mol().v() -= 2*vn*normal_;
                    }
                }
            }
        }
    }
}

void polyInternalWall::controlAfterVelocityII()
{}

void polyInternalWall::output
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


void polyInternalWall::calculateProperties()
{}

void polyInternalWall::updateProperties(const dictionary& newDict)
{
    //- the main controller properties should be updated first
    updateStateControllerProperties(newDict);

    propsDict_ = newDict.subDict(typeName + "Properties");


    readProperties();
}


void polyInternalWall::readProperties()
{

}



} // End namespace Foam

// ************************************************************************* //
