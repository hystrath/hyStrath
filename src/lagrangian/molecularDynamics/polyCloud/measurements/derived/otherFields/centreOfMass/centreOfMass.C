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

#include "centreOfMass.H"
#include "addToRunTimeSelectionTable.H"
#include "fvc.H"
#include "IFstream.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(centreOfMass, 0);

addToRunTimeSelectionTable(polyField, centreOfMass, dictionary);



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
centreOfMass::centreOfMass
(
    Time& t,
    const polyMesh& mesh,
    polyMoleculeCloud& molCloud,
    const dictionary& dict
)
:
    polyField(t, mesh, molCloud, dict),
    propsDict_(dict.subDict(typeName + "Properties")),
    fieldName_(propsDict_.lookup("fieldName")),
    fields_(t, mesh, "dummy"),
    molIds_()
{
   // choose molecule ids to sample

    molIds_.clear();

    selectIds ids
    (
        molCloud_.cP(),
        propsDict_
    );

    molIds_ = ids.molIds();

}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

centreOfMass::~centreOfMass()
{}





// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void centreOfMass::createField()
{
    calculateField();
}


//- call this function every time-step before the state and flux objects are cleaned
void centreOfMass::calculateField()
{

    IDLList<polyMolecule>::iterator mol(molCloud_.begin());
    
    vector centre = vector::zero;
    scalar nMols = 0.0;
    
    for (mol = molCloud_.begin(); mol != molCloud_.end(); ++mol)
    {
        if(findIndex(molIds_, mol().id()) != -1)
        {
            centre += mol().position();
            nMols += 1.0;
        }
    }    
    
    if(Pstream::parRun())
    {
        reduce(centre, sumOp<vector>());
        reduce(nMols, sumOp<scalar>());
    }
    
    if(nMols > 0)
    {
        centre /= nMols;
    }

    if(Pstream::master())
    {    
        centreOfMass_.append(centre);
    }
}
    
void centreOfMass::writeField()
{
    const Time& runTime = time_.time();

    if(runTime.outputTime())
    {
        if(Pstream::master())
        {
            vectorField centreOfMass (centreOfMass_.size(), vector::zero);
            scalarField timeField (centreOfMass_.size(), 0.0);
            
            centreOfMass.transfer(centreOfMass_);
            centreOfMass_.clear();

            const scalar& deltaT = time_.time().deltaT().value();
            
            forAll(timeField, i)
            {
                timeField[timeField.size()-i-1]=runTime.timeOutputValue()-(deltaT*i);
            }            
            
            writeTimeData
            (
                casePath_,
                "centreOfMass_"+fieldName_+".xyz",
                timeField,
                centreOfMass,
                true
            );
        }
        
        centreOfMass_.clear();        
    }
}

void centreOfMass::measureDuringForceComputation
(
    polyMolecule* molI,
    polyMolecule* molJ
){}

void centreOfMass::measureDuringForceComputationSite
(   
    polyMolecule* molI,
    polyMolecule* molJ,
    label sI,
    label sJ
){}

const propertyField& centreOfMass::fields() const
{
    return  fields_;
}



} // End namespace Foam

// ************************************************************************* //
