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

#include "polyForceTwoShear.H"
#include "addToRunTimeSelectionTable.H"
#include "IFstream.H"
#include "graph.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(polyForceTwoShear, 0);
addToRunTimeSelectionTable(polyStateController, polyForceTwoShear, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
polyForceTwoShear::polyForceTwoShear
(
    Time& t,
    polyMoleculeCloud& molCloud,
    const dictionary& dict
)
:
    polyStateController(t, molCloud, dict),
    propsDict_(dict.subDict(typeName + "Properties")),

    molIds_(),
    forceMagLiquid_(propsDict_.lookup("forceLiquid"))
{

    writeInTimeDir_ = true;
    writeInCase_ = true;

//     singleValueController() = true;
    
    setBoundBox(propsDict_, bb_liquid_A_, "liquidA");
    setBoundBox(propsDict_, bb_liquid_B_, "liquidB");    
    setBoundBox(propsDict_, bb_gas_A_, "gasA");
    setBoundBox(propsDict_, bb_gas_B_, "gasB");

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

polyForceTwoShear::~polyForceTwoShear()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void polyForceTwoShear::initialConfiguration()
{}

void polyForceTwoShear::controlBeforeVelocityI()
{}


void polyForceTwoShear::controlBeforeMove()
{}

void polyForceTwoShear::controlBeforeForces()
{}

void polyForceTwoShear::controlDuringForces
(
    polyMolecule* molI,
    polyMolecule* molJ
)
{}

void polyForceTwoShear::controlAfterForces()
{
    scalar molsLiquid = 0.0;
    scalar molsGas = 0.0;     
    
    // measurements of instantaneous density
    {
    
        IDLList<polyMolecule>::iterator mol(molCloud_.begin());

        for (mol = molCloud_.begin(); mol != molCloud_.end(); ++mol)
        {
            if(bb_liquid_A_.contains(mol().position()) || bb_liquid_B_.contains(mol().position()))
            {
                if(findIndex(molIds_, mol().id()) != -1)
                {
                    molsLiquid += 1.0;
                }
            }
            
            if(bb_gas_A_.contains(mol().position()) || bb_gas_B_.contains(mol().position()))
            {
                if(findIndex(molIds_, mol().id()) != -1)
                {
                    molsGas += 1.0;
                }
            }            
        }
        
        // parallel processing
        
        if(Pstream::parRun())
        {
            reduce(molsLiquid, sumOp<scalar>());
            reduce(molsGas, sumOp<scalar>());
        }
    }

    forceMagGas_ = forceMagLiquid_*molsLiquid/molsGas;
    
    Info << "polyForceTwoShear: control - force liquid = "  << forceMagLiquid_ 
         << " , force gas = " << forceMagGas_
         << endl;

    IDLList<polyMolecule>::iterator mol(molCloud_.begin());

    vector totalForceLiquid = vector::zero;
    vector totalForceGas = vector::zero;
    
    for (mol = molCloud_.begin(); mol != molCloud_.end(); ++mol)
    {
        if(bb_liquid_A_.contains(mol().position()) || bb_liquid_B_.contains(mol().position()))
        {
            if(findIndex(molIds_, mol().id()) != -1)
            {
                const scalar& massI = molCloud_.cP().mass(mol().id());

                vector force = forceMagLiquid_;
                
                mol().a() += force/massI;

                totalForceLiquid += force;
            }
        }
        
        if(bb_gas_A_.contains(mol().position()) || bb_gas_B_.contains(mol().position()))
        {
            if(findIndex(molIds_, mol().id()) != -1)
            {
                const scalar& massI = molCloud_.cP().mass(mol().id());
                
                vector force = -forceMagGas_;
                
                mol().a() += force/massI;
                
                totalForceGas += force;
            }
        }        
        
    }
    
    if(Pstream::parRun())
    {
        reduce(totalForceLiquid, sumOp<vector>());
        reduce(totalForceGas, sumOp<vector>());
    }    
	
    
    Info << "polyForceTwoShear: Conservation- total force liquid = "  << totalForceLiquid 
         << " , total force gas = " << totalForceGas
         << " , diff = " << totalForceLiquid + totalForceGas 
         << endl;
}

void polyForceTwoShear::controlAfterVelocityII()
{}

void polyForceTwoShear::calculateProperties()
{}

void polyForceTwoShear::output
(
    const fileName& fixedPathName,
    const fileName& timePath
)
{}

void polyForceTwoShear::updateProperties(const dictionary& newDict)
{
    //- the main controller properties should be updated first
    updateStateControllerProperties(newDict);

    propsDict_ = newDict.subDict(typeName + "Properties");

    readProperties();
}

void polyForceTwoShear::readProperties()
{}

void polyForceTwoShear::setBoundBox
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

} // End namespace Foam

// ************************************************************************* //
