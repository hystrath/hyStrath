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

#include "fadeInitialise.H"
#include "addToRunTimeSelectionTable.H"
#include "IFstream.H"
#include "graph.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(fadeInitialise, 0);

addToRunTimeSelectionTable(polyStateController, fadeInitialise, dictionary);





// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
fadeInitialise::fadeInitialise
(
    Time& t,
    polyMoleculeCloud& molCloud,
    const dictionary& dict
)
:
    polyStateController(t,  molCloud, dict),
    propsDict_(dict.subDict(typeName + "Properties")),
    molIds_()
{
    writeInTimeDir_ = false;
    writeInCase_ = false;


    molIds_.clear();

    selectIds ids
    (
        molCloud_.cP(),
        propsDict_
    );

    molIds_ = ids.molIds();

    n_ = readLabel(propsDict_.lookup("n"));
    tauT_ = readScalar(propsDict_.lookup("tauT"));   
    deltaT_ = time_.time().deltaT().value();
    t_ = 0.0;
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

fadeInitialise::~fadeInitialise()
{}



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void fadeInitialise::initialConfiguration()
{
    IDLList<polyMolecule>::iterator mol(molCloud_.begin());

    for (mol = molCloud_.begin(); mol != molCloud_.end(); ++mol)
    {
        if(findIndex(molIds_, mol().id()) != -1)
        {         
            mol().fraction() = 0.0;
        }
    }
}

void fadeInitialise::controlBeforeVelocityI()
{}

void fadeInitialise::controlBeforeMove()
{}

void fadeInitialise::controlBeforeForces()
{}

void fadeInitialise::controlDuringForces
(
    polyMolecule* molI,
    polyMolecule* molJ
)
{}

void fadeInitialise::controlAfterForces()
{

}


void fadeInitialise::controlAfterVelocityII()
{}

void fadeInitialise::calculateProperties()
{
    t_ += deltaT_;
    
    if(t_ < tauT_)
    {
        IDLList<polyMolecule>::iterator mol(molCloud_.begin());

        for (mol = molCloud_.begin(); mol != molCloud_.end(); ++mol)
        {
            if(findIndex(molIds_, mol().id()) != -1)
            {    
                if(t_ < 0.5*tauT_)
                {
                    mol().fraction() = 0.5*Foam::pow((2.0*t_/tauT_), scalar(n_));
                }
                else
                {
                    mol().fraction() = 1.0 
                            - 0.5*mag(Foam::pow((2.0*(t_-tauT_)/tauT_), scalar(n_)));
                }
            }
        }
    }
    else
    {
        Info << "End of equilibration" << endl;

        IDLList<polyMolecule>::iterator mol(molCloud_.begin());
        
        for (mol = molCloud_.begin(); mol != molCloud_.end(); ++mol)
        {
            if(findIndex(molIds_, mol().id()) != -1)
            {        
                mol().fraction() = 1.0;
            }
        }
    }

}

void fadeInitialise::output
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

void fadeInitialise::updateProperties(const dictionary& newDict)
{
    //- the main controller properties should be updated first
    updateStateControllerProperties(newDict);

    propsDict_ = newDict.subDict(typeName + "Properties");

}






} // End namespace Foam

// ************************************************************************* //
