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

#include "polyARCHER.H"
#include "addToRunTimeSelectionTable.H"
#include "graph.H"
#include "OFstream.H"
#include "IFstream.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(polyARCHER, 0);

addToRunTimeSelectionTable(polyField, polyARCHER, dictionary);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
polyARCHER::polyARCHER
(
    Time& t,
    const polyMesh& mesh,
    polyMoleculeCloud& molCloud,
    const dictionary& dict
)
:
    polyField(t, mesh, molCloud, dict),
    propsDict_(dict.subDict(typeName + "Properties")),    
    fields_(t, mesh, "dummy"),
    nameOfFile_(propsDict_.lookup("fileName"))    
{
    pathName_ = time_.time().path();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

polyARCHER::~polyARCHER()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void polyARCHER::createField()
{
    outputInitialisation();
}

scalar polyARCHER::getTotalEnergy()
{
    scalar kE = 0.0;
    scalar angularKE = 0.0;
    scalar PE = 0.0;

    label nMols = molCloud_.size();

    {
        IDLList<polyMolecule>::iterator mol(molCloud_.begin());
    
        for
        (
            mol = molCloud_.begin();
            mol != molCloud_.end();
            ++mol
        )
        {
            label molId = mol().id();
    
            scalar molMass(molCloud_.cP().mass(molId));

            const vector& molV(mol().v());

            const diagTensor& molMoI(molCloud_.cP().momentOfInertia(mol().id()));

            const vector& molOmega(inv(molMoI) & mol().pi());
    
            kE += 0.5*molMass*magSqr(molV);
            angularKE += 0.5*(molOmega & molMoI & molOmega);
            PE += mol().potentialEnergy();
        }
    }

    if (Pstream::parRun())
    {
        reduce(kE, sumOp<scalar>());
        reduce(angularKE, sumOp<scalar>());
        reduce(PE, sumOp<scalar>());
        reduce(nMols, sumOp<label>());
    }

    scalar totalEnergy = 0.0;
    
    if(nMols > 0)
    {
        totalEnergy = (kE + angularKE + PE)/nMols;
    }            
    
    return totalEnergy;
}

void polyARCHER::calculateField()
{
    if(time_.outputTime())
    {
        outputTime();
    }
}

void polyARCHER::writeField()
{

}

void polyARCHER::measureDuringForceComputation
(
    polyMolecule* molI,
    polyMolecule* molJ
)
{}

void polyARCHER::measureDuringForceComputationSite
(
    polyMolecule* molI,
    polyMolecule* molJ,
    label sI,
    label sJ
)
{}

const propertyField& polyARCHER::fields() const
{
    return  fields_;
}

} // End namespace Foam

// ************************************************************************* //
