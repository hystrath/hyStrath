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

#include "polyMoleculePotentialEnergy.H"
#include "addToRunTimeSelectionTable.H"
#include "graph.H"
// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(polyMoleculePotentialEnergy, 0);
addToRunTimeSelectionTable(polyField, polyMoleculePotentialEnergy, dictionary);

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
polyMoleculePotentialEnergy::polyMoleculePotentialEnergy
(
    Time& t,
    const polyMesh& mesh,
    polyMoleculeCloud& molCloud,
    const dictionary& dict
)
:
    polyField(t, mesh, molCloud, dict),
    fields_(t, mesh, "dummy"),
    propsDict_(dict.subDict(typeName + "Properties")),
    molPoint_(propsDict_.lookup("molPoint")),
    cellI_(-1),
    molTrackingNumber_(-1)
{
    cellI_ = mesh_.findCell(molPoint_);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

polyMoleculePotentialEnergy::~polyMoleculePotentialEnergy()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void polyMoleculePotentialEnergy::createField()
{
    if(cellI_ != -1)
    {
        const List<polyMolecule*>& molsInCell = molCloud_.cellOccupancy()[cellI_];

        scalar rD = GREAT;

        forAll(molsInCell, m)
        {
            polyMolecule* molI = molsInCell[m];

            scalar rIPMag = mag(molI->position() - molPoint_);

            if(rIPMag < rD)
            {
                rD = rIPMag;
                molTrackingNumber_ = molI->trackingNumber();
            }
        }

        Pout << "tracking number: " << molTrackingNumber_  << endl;

        forAll(molsInCell, m)
        {
            polyMolecule* molI = molsInCell[m];

            if(molTrackingNumber_ == molI->trackingNumber())
            {
                if(Pstream::parRun())
                {     
                    Pout << "molecule at position: " << molI->position() 
                        << ", PE: " << molI->potentialEnergy()
                        << endl;
                }
                else
                {
                    Info << "molecule at position: " << molI->position() 
                        << ", PE: " << molI->potentialEnergy()
                        << endl;
                }
            }
        }
    }
}


void polyMoleculePotentialEnergy::calculateField()
{
    IDLList<polyMolecule>::iterator mol(molCloud_.begin());

    bool notFound = true;

    for (mol = molCloud_.begin(); mol != molCloud_.end(); ++mol)
    {
        if(mol().trackingNumber() == molTrackingNumber_)
        {
            cellI_ = mol().cell();

            notFound = false;

            if(Pstream::parRun())
            {     
                Pout << "molecule at position: " << mol().position() 
                    << ", PE: " << mol().potentialEnergy()
                    << ", current processor: " << Pstream::myProcNo()
                    << endl;
            }
            else
            {
                Info << "molecule at position: " << mol().position() 
                    << ", PE: " << mol().potentialEnergy()
                    << endl;
            }
        }
    }

    if(notFound)
    {
        cellI_ = -1;
    }
}

void polyMoleculePotentialEnergy::writeField()
{}

void polyMoleculePotentialEnergy::measureDuringForceComputation
(
    polyMolecule* molI,
    polyMolecule* molJ
)
{}

void polyMoleculePotentialEnergy::measureDuringForceComputationSite
(
    polyMolecule* molI,
    polyMolecule* molJ,
    label sI,
    label sJ
)
{}

const propertyField& polyMoleculePotentialEnergy::fields() const
{
    return  fields_;
}

} // End namespace Foam

// ************************************************************************* //
