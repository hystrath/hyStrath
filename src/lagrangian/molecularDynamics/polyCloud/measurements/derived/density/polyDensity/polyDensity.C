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

#include "polyDensity.H"
#include "addToRunTimeSelectionTable.H"
#include "fvc.H"
#include "IFstream.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(polyDensity, 0);

addToRunTimeSelectionTable(polyField, polyDensity, dictionary);



void polyDensity::setBoundBox
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


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
polyDensity::polyDensity
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
    
    setBoundBox(propsDict_, bb_, "samplingRegion");  
    
    volume_ = bb_.span().x() * bb_.span().y() * bb_.span().z();    
    
    Info << "Volume = " << volume_ << endl;
    
    
    if (propsDict_.found("volume"))
    {
        scalar volume  = readScalar(propsDict_.lookup("volume"));
        
        volume_ = volume;
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

polyDensity::~polyDensity()
{}





// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void polyDensity::createField()
{
}
//- call this function every time-step before the state and flux objects are cleaned
void polyDensity::calculateField()
{
    IDLList<polyMolecule>::iterator mol(molCloud_.begin());
    
    scalar mass = 0.0;
    
    for (mol = molCloud_.begin(); mol != molCloud_.end(); ++mol)
    {
        if(bb_.contains(mol().position()))
        {
            if(findIndex(molIds_, mol().id()) != -1)
            {
                const scalar& massI = molCloud_.cP().mass(mol().id());
                mass += massI;            
            }
        }
    }

    if(Pstream::parRun())
    {
        reduce(mass, sumOp<scalar>());
    }
    
    scalar density = mass/volume_;
    
    densities_.append(density);

}

void polyDensity::writeField()
{
    const Time& runTime = time_;

    if(runTime.outputTime())
    {
        if(Pstream::master())
        {
            densities_.shrink();

            scalarField timeField (densities_.size(), 0.0);
            scalarField densities (densities_.size(), 0.0);
            
            densities.transfer(densities_);
            densities_.clear();

            const scalar& deltaT = time_.deltaT().value();
            
            forAll(timeField, i)
            {
                timeField[timeField.size()-i-1]=runTime.timeOutputValue()-(deltaT*i);
            }
            
            writeTimeData
            (
                casePath_,
                "density_bb_"+fieldName_+".xy",
                timeField,
                densities,
                true
            );

            const reducedUnits& rU = molCloud_.redUnits();
    
            if(rU.outputSIUnits())
            {
                writeTimeData
                (
                    casePath_,
                    "density_bb_"+fieldName_+"_SI.xy",
                    timeField*rU.refTime(),
                    densities*rU.refMassDensity(),
                    true
                );
            }
        }
    }
}

void polyDensity::measureDuringForceComputation
(
    polyMolecule* molI,
    polyMolecule* molJ
){}

void polyDensity::measureDuringForceComputationSite
(   
    polyMolecule* molI,
    polyMolecule* molJ,
    label sI,
    label sJ
){}

const propertyField& polyDensity::fields() const
{
    return  fields_;
}

}// End namespace Foam

// ************************************************************************* //
