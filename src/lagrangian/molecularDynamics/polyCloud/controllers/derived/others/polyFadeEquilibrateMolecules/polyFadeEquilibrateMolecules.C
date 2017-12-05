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

#include "polyFadeEquilibrateMolecules.H"
#include "addToRunTimeSelectionTable.H"
#include "IFstream.H"
#include "graph.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(polyFadeEquilibrateMolecules, 0);
addToRunTimeSelectionTable(polyStateController, polyFadeEquilibrateMolecules, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
polyFadeEquilibrateMolecules::polyFadeEquilibrateMolecules
(
    Time& t,
    polyMoleculeCloud& molCloud,
    const dictionary& dict
)
:
    polyStateController(t, molCloud, dict),
    propsDict_(dict.subDict(typeName + "Properties")),

    molIds_()
//     nTimeSteps_(0.0),
//     force_(vector::zero)
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
    
    tauT_ = readScalar(propsDict_.lookup("tauT"));

    if(tauT_ > 0.0)
    {
        n_ = readLabel(propsDict_.lookup("n"));
    }
    
    deltaT_ = time_.deltaT().value(); 
    t_ = 0.0;
    
    readProperties();
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

polyFadeEquilibrateMolecules::~polyFadeEquilibrateMolecules()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void polyFadeEquilibrateMolecules::initialConfiguration()
{
    IDLList<polyMolecule>::iterator mol(molCloud_.begin());

    for (mol = molCloud_.begin(); mol != molCloud_.end(); ++mol)
    {
        if(findIndex(molIds_, mol().id()) != -1)
        {
            mol().fraction() = 0;
        }
    }
}

void polyFadeEquilibrateMolecules::controlBeforeVelocityI()
{}

void polyFadeEquilibrateMolecules::controlBeforeMove()
{}

void polyFadeEquilibrateMolecules::controlBeforeForces()
{}

void polyFadeEquilibrateMolecules::controlDuringForces
(
    polyMolecule* molI,
    polyMolecule* molJ
)
{}

void polyFadeEquilibrateMolecules::controlAfterForces()
{
    IDLList<polyMolecule>::iterator mol(molCloud_.begin());
    
    t_ += deltaT_;
    
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

void polyFadeEquilibrateMolecules::controlAfterVelocityII()
{}

void polyFadeEquilibrateMolecules::calculateProperties()
{}

void polyFadeEquilibrateMolecules::output
(
    const fileName& fixedPathName,
    const fileName& timePath
)
{

}

void polyFadeEquilibrateMolecules::updateProperties(const dictionary& newDict)
{
    //- the main controller properties should be updated first
    updateStateControllerProperties(newDict);

    propsDict_ = newDict.subDict(typeName + "Properties");


    readProperties();
}

void polyFadeEquilibrateMolecules::readProperties()
{}

} // End namespace Foam

// ************************************************************************* //
