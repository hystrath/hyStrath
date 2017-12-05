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

#include "polyBoundedSwitchIds.H"
#include "addToRunTimeSelectionTable.H"
#include "IFstream.H"
#include "graph.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(polyBoundedSwitchIds, 0);

addToRunTimeSelectionTable(polyStateController, polyBoundedSwitchIds, dictionary);




// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
polyBoundedSwitchIds::polyBoundedSwitchIds
(
    Time& t,
    polyMoleculeCloud& molCloud,
    const dictionary& dict
)
:
    polyStateController(t,  molCloud, dict),
    propsDict_(dict.subDict(typeName + "Properties"))   
{
    writeInTimeDir_ = false;
    writeInCase_ = false;

    {
        initialMolIds_.clear();

        selectIds ids
        (
            molCloud_.cP(),
            propsDict_,
            "initialIds"
        );

        initialMolIds_ = ids.molIds();
    }
    
    {
        finalMolIds_.clear();

        selectIds ids
        (
            molCloud_.cP(),
            propsDict_,
            "finalIds"
        );

        finalMolIds_ = ids.molIds();
    }    
    
    readProperties();

    vector startPoint = propsDict_.lookup("startPoint");
    vector endPoint = propsDict_.lookup("endPoint");
    box_.resetBoundedBox(startPoint, endPoint);
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

polyBoundedSwitchIds::~polyBoundedSwitchIds()
{}



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void polyBoundedSwitchIds::initialConfiguration()
{}

void polyBoundedSwitchIds::controlBeforeVelocityI()
{}

void polyBoundedSwitchIds::controlBeforeMove()
{}

void polyBoundedSwitchIds::controlBeforeForces()
{}

void polyBoundedSwitchIds::controlDuringForces
(
    polyMolecule* molI,
    polyMolecule* molJ
)
{}

void polyBoundedSwitchIds::controlAfterForces()
{
    Info << "polyBoundedSwitchIds: control" << endl;

    IDLList<polyMolecule>::iterator mol(molCloud_.begin());
    
    label nMols = 0;
    
    for (mol = molCloud_.begin(); mol != molCloud_.end(); ++mol)
    {
        if(findIndex(initialMolIds_, mol().id()) != -1)
        {
            if(box_.contains(mol().position()))
            {
                mol().id() = finalMolIds_[0];
                nMols++;
            }
        }
    }
    

    if(Pstream::parRun())
    {
        reduce(nMols, sumOp<label>());
    }    
    
    Info << "number of molecules switched = " << nMols << endl;
}


void polyBoundedSwitchIds::controlAfterVelocityII()
{}

void polyBoundedSwitchIds::calculateProperties()
{}

void polyBoundedSwitchIds::output
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

void polyBoundedSwitchIds::updateProperties(const dictionary& newDict)
{
    //- the main controller properties should be updated first
    updateStateControllerProperties(newDict);

    propsDict_ = newDict.subDict(typeName + "Properties");

//     model_->updateProperties(propsDict_);

//     readProperties();
}


void polyBoundedSwitchIds::readProperties()
{

}



} // End namespace Foam

// ************************************************************************* //
