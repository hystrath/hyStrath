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

#include "pressureFlux.H"
#include "addToRunTimeSelectionTable.H"
#include "IFstream.H"
#include "graph.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(pressureFlux, 0);

addToRunTimeSelectionTable(polyStateController, pressureFlux, dictionary);



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

void pressureFlux::setBoundBox
(
    const dictionary& propsDict,
    boundedBox& bb,
    const word& name 
)
{
    const dictionary& dict(propsDict.subDict(name));
    
    vector startPoint = dict.lookup("startPoint");
    vector endPoint = dict.lookup("endPoint");

    bb.resetBoundedBox(startPoint, endPoint);
}

// Construct from components
pressureFlux::pressureFlux
(
    Time& t,
    polyMoleculeCloud& molCloud,
    const dictionary& dict
)
:
    polyStateController(t,  molCloud, dict),
    propsDict_(dict.subDict(typeName + "Properties")),
    fileName_(propsDict_.lookup("fileName")),  
    molIds_()
//     nTimeSteps_(0.0),
//     force_(vector::zero)
   
{
    writeInTimeDir_ = true;
    writeInCase_ = true;

    
    molIds_.clear();

    selectIds ids
    (
        molCloud_.cP(),
        propsDict_
    );

    molIds_ = ids.molIds();

    
    setBoundBox(propsDict_, box_, "controlBoundBox");
    
    d_ = propsDict_.lookup("forceDirection");
    d_ /= mag(d_);
    
    area_ = readScalar(propsDict_.lookup("area"));

    targetPressure_ = readScalar(propsDict_.lookup("pressureMPa"));    
    
    targetPressure_ *= 1e6;
    
    targetPressure_ /= molCloud_.redUnits().refPressure();
    
    Info << "target Pressure RU " <<  targetPressure_ << endl;

    min_ = 100.0; // default
    
    if (propsDict_.found("minMols"))
    {
        min_ = readScalar(propsDict_.lookup("minMols"));         
    }
    
    force_ = vector::zero;
    nTimeSteps_ = 0.0;
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

pressureFlux::~pressureFlux()
{}



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void pressureFlux::initialConfiguration()
{}

void pressureFlux::controlBeforeVelocityI()
{}

void pressureFlux::controlBeforeMove()
{}

void pressureFlux::controlBeforeForces()
{}

void pressureFlux::controlDuringForces
(
    polyMolecule* molI,
    polyMolecule* molJ
)
{}

void pressureFlux::controlAfterForces()
{
    Info << "pressureFlux: control" << endl;

    // measure number of molecules in box 
    
    scalar nMols = 0.0;
    {
        IDLList<polyMolecule>::iterator mol(molCloud_.begin());
        
        for (mol = molCloud_.begin(); mol != molCloud_.end(); ++mol)
        {
            if(findIndex(molIds_, mol().id()) != -1)
            {
                if(box_.contains(mol().position()))
                {
                    nMols += 1.0;
                }
            }
        }    
    }
    
    if(Pstream::parRun())
    {
        reduce(nMols, sumOp<scalar>());
    }
    
    nMols_ += nMols;
    
    vector force = vector::zero;
    
    if(nMols > min_)
    {
        force = targetPressure_*area_*d_/nMols;
    }
    else
    {
        FatalErrorIn("polyBinsMethod::polyBinsMethod()")
            << "Too few number of molecules inside control region = "
            << nMols 
            << nl 
            << exit(FatalError);         
    }

    IDLList<polyMolecule>::iterator mol(molCloud_.begin());

    for (mol = molCloud_.begin(); mol != molCloud_.end(); ++mol)
    {
        if(findIndex(molIds_, mol().id()) != -1)
        {
            if(box_.contains(mol().position()))
            {
                const scalar& massI = molCloud_.cP().mass(mol().id());
                
                mol().a() += force/massI;

                force_ += force;
            }
        }
    }
    
    nTimeSteps_ += 1.0;
    
}


void pressureFlux::controlAfterVelocityII()
{}

void pressureFlux::calculateProperties()
{}

void pressureFlux::output
(
    const fileName& fixedPathName,
    const fileName& timePath
)
{
    const Time& runTime = time_.time();

    if(runTime.outputTime())
    {
        if(Pstream::master())
        {
            scalarField timeField(1, time_.time().timeOutputValue());
            vectorField force(1, force_/nTimeSteps_);            
            scalarField nMols(1, nMols_/nTimeSteps_);
            
            writeTimeData
            (
                fixedPathName,
                "pressureFlux_"+fileName_+"_force.xyz",
                timeField,
                force,
                true
            );
            
            writeTimeData
            (
                fixedPathName,
                "pressureFlux_"+fileName_+"_N.xyz",
                timeField,
                nMols,
                true
            );            
        }
        
        force_ = vector::zero;
        nTimeSteps_ = 0.0;
    }
    
    
}

void pressureFlux::updateProperties(const dictionary& newDict)
{
    //- the main controller properties should be updated first
    updateStateControllerProperties(newDict);

    propsDict_ = newDict.subDict(typeName + "Properties");

//     model_->updateProperties(propsDict_);

//     readProperties();
}




} // End namespace Foam

// ************************************************************************* //
