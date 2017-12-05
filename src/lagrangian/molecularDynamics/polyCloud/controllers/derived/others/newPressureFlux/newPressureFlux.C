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

#include "newPressureFlux.H"
#include "addToRunTimeSelectionTable.H"
#include "IFstream.H"
#include "graph.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(newPressureFlux, 0);

addToRunTimeSelectionTable(polyStateController, newPressureFlux, dictionary);



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

void newPressureFlux::setBoundBox
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
newPressureFlux::newPressureFlux
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
    
    deltaY_ = 0.0;
    
    if (propsDict_.found("deltaY"))
    {
        deltaY_ = readScalar(propsDict_.lookup("deltaY"));         
    }    
    
//     force_ = vector::zero;
//     nTimeSteps_ = 0.0;
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

newPressureFlux::~newPressureFlux()
{}



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void newPressureFlux::initialConfiguration()
{}

void newPressureFlux::controlBeforeVelocityI()
{}

void newPressureFlux::controlBeforeMove()
{}

void newPressureFlux::controlBeforeForces()
{}

void newPressureFlux::controlDuringForces
(
    polyMolecule* molI,
    polyMolecule* molJ
)
{}

void newPressureFlux::controlAfterForces()
{
    Info << "newPressureFlux: control" << endl;

    // measure number of molecules in the entire water slab 
    
    vector centre = vector::zero;
    
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
                    centre += mol().position();
                }
            }
        }    
    }
    
    if(Pstream::parRun())
    {
        reduce(nMols, sumOp<scalar>());
        reduce(centre, sumOp<vector>());
    }
    
    if(nMols > 0.0)
    {
        centre /= nMols;
        centre_.append(centre);
    }

    
    // measure molecules above the centre point of water slab 
    scalar mols = 0.0;
    
    {
        IDLList<polyMolecule>::iterator mol(molCloud_.begin());
        
        for (mol = molCloud_.begin(); mol != molCloud_.end(); ++mol)
        {
            if(findIndex(molIds_, mol().id()) != -1)
            {
                if(box_.contains(mol().position()))
                {
                    if(mol().position().y() > (centre.y() + deltaY_))
                    {
                        mols += 1.0;
                    }
                }
            }
        }    
    }    
    
    if(Pstream::parRun())
    {
        reduce(mols, sumOp<scalar>());
    }
    
    vector force = vector::zero;
    
    if(mols > min_)
    {
        force = targetPressure_*area_*d_/mols;
    }
    else
    {
        FatalErrorIn("polyBinsMethod::polyBinsMethod()")
            << "Too few number of molecules inside control region = "
            << mols 
            << nl 
            << exit(FatalError);         
    }
    
    Info << "FORCE = " << force << endl;
    
    force_.append(force);
    mols_.append(mols);
    
    vector totalForce = vector::zero;
    
    // apply force to those molecules above the centre point of water slab 
    {
        IDLList<polyMolecule>::iterator mol(molCloud_.begin());

        for (mol = molCloud_.begin(); mol != molCloud_.end(); ++mol)
        {
            if(findIndex(molIds_, mol().id()) != -1)
            {
                if(box_.contains(mol().position()))
                {
                    if(mol().position().y() > (centre.y() + deltaY_))
                    {                
                        const scalar& massI = molCloud_.cP().mass(mol().id());
                    
                        mol().a() += force/massI;
                        
                        totalForce += force;
                    }
                }
            }
        }
    }
    
    
    if(Pstream::parRun())
    {
        reduce(totalForce, sumOp<vector>());
    }
    
    totalForce_.append(totalForce);
}


void newPressureFlux::controlAfterVelocityII()
{}

void newPressureFlux::calculateProperties()
{}

void newPressureFlux::output
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
            scalarField mols (mols_.size(), 0.0);
            vectorField force (force_.size(), vector::zero);
            vectorField totalForce (totalForce_.size(), vector::zero);                        
            vectorField centre (centre_.size(), vector::zero);
            scalarField timeField (centre_.size(), 0.0);
            
            
            mols.transfer(mols_);            
            force.transfer(force_);            
            totalForce.transfer(totalForce_);             
            
            centre.transfer(centre_);
            

            const scalar& deltaT = time_.time().deltaT().value();
            
            forAll(timeField, i)
            {
                timeField[timeField.size()-i-1]=runTime.timeOutputValue()-(deltaT*i);
            }             
            
            writeTimeData
            (
                fixedPathName,
                "newPressureFlux_"+fileName_+"_centre.xyz",
                timeField,
                centre,
                true
            );
            
            writeTimeData
            (
                fixedPathName,
                "newPressureFlux_"+fileName_+"_forcePerMol.xyz",
                timeField,
                force,
                true
            );

            
            writeTimeData
            (
                fixedPathName,
                "newPressureFlux_"+fileName_+"_totalForce.xyz",
                timeField,
                totalForce,
                true
            );
            
            writeTimeData
            (
                fixedPathName,
                "newPressureFlux_"+fileName_+"_N.xyz",
                timeField,
                mols,
                true
            );            
        }
        
        mols_.clear();            
        force_.clear();            
        totalForce_.clear();           
        centre_.clear();        
    }
}

void newPressureFlux::updateProperties(const dictionary& newDict)
{
    //- the main controller properties should be updated first
    updateStateControllerProperties(newDict);

    propsDict_ = newDict.subDict(typeName + "Properties");

//     model_->updateProperties(propsDict_);

//     readProperties();
}




} // End namespace Foam

// ************************************************************************* //
