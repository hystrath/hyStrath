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

#include "polyVelocityBins.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(polyVelocityBins, 0);

addToRunTimeSelectionTable(polyStateController, polyVelocityBins, dictionary);



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
polyVelocityBins::polyVelocityBins
(
    Time& t,
    polyMoleculeCloud& molCloud,
    const dictionary& dict
)
:
    polyStateController(t, molCloud, dict),
    propsDict_(dict.subDict(typeName + "Properties")),
    lambda_(readScalar(propsDict_.lookup("lambda"))),
    componentControl_(false),
    X_(false),
    Y_(false),
    Z_(false),
    molIds_()
{
    writeInTimeDir_ = true;
    writeInCase_ = true;

//     singleValueController() = true;

    velocity_ = propsDict_.lookup("velocity");

    if (propsDict_.found("componentControl"))
    {
        componentControl_ = Switch(propsDict_.lookup("componentControl"));

        if(componentControl_)
        {
            if (propsDict_.found("X"))
            {
                X_ = Switch(propsDict_.lookup("X"));
            }
    
            if (propsDict_.found("Y"))
            {
                Y_ = Switch(propsDict_.lookup("Y"));
            }
    
            if (propsDict_.found("Z"))
            {
                Z_ = Switch(propsDict_.lookup("Z"));
            }

            if(!X_ && !Y_ && !Z_)
            {
                FatalErrorIn("polyVelocityBins::polyVelocityBins()")
                    << "At least one component (X, Y, Z) should be chosen " << nl << "in: "
                    << time_.time().system()/"controllersDict"
                    << exit(FatalError);
            }
        }
    }
    

    molIds_.clear();

    selectIds ids
    (
        molCloud_.cP(),
        propsDict_
    );

    molIds_ = ids.molIds();

    binModel_ =  autoPtr<binModel>
    (
        binModel::New(mesh_, propsDict_)
    );

}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

polyVelocityBins::~polyVelocityBins()
{}



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void polyVelocityBins::initialConfiguration()
{
}

void polyVelocityBins::controlBeforeVelocityI()
{
    const List< DynamicList<polyMolecule*> >& cellOccupancy
        = molCloud_.cellOccupancy();

    const labelList& cells = mesh_.cellZones()[regionId_];
    
    const label& nBins = binModel_->nBins();
    
    scalarField mass(nBins, 0.0);
    vectorField momentum(nBins, vector::zero);

    forAll(cells, c)
    {
        const label& cell = cells[c];
        const List<polyMolecule*>& molsInCell = cellOccupancy[cell];

        forAll(molsInCell, mIC)
        {
            polyMolecule* molI = molsInCell[mIC];

            const vector& rI = molI->position();

            label n = binModel_->isPointWithinBin(rI, cell);

            if(n != -1)
            {
                if(findIndex(molIds_, molI->id()) != -1)
                {
                    const scalar& massI = molCloud_.cP().mass(molI->id());

                    mass[n] += massI;
                    momentum[n] += massI*molI->v();
                }
            }
        }
    }

    if(Pstream::parRun())
    {
        forAll(mass, i)
        {
            reduce(momentum[i], sumOp<vector>());
            reduce(mass[i], sumOp<scalar>());
        }
    }

    vectorField velMeasured(nBins, vector::zero);
    vectorField deltaU(nBins, vector::zero);
    
    forAll(velMeasured, n)
    {
        if(mass[n] > 0.0)
        {
            velMeasured[n] = momentum[n]/mass[n];
            deltaU[n] = (velocity_ - velMeasured[n])*lambda_;    
//             Info << "vel = " << velMeasured[n] << endl;
        }
    }
    
    // control 
    
    forAll(cells, c)
    {
        const label& cell = cells[c];
        const List<polyMolecule*>& molsInCell = cellOccupancy[cell];

        forAll(molsInCell, mIC)
        {
            polyMolecule* molI = molsInCell[mIC];

            const vector& rI = molI->position();

            label n = binModel_->isPointWithinBin(rI, cell);

            if(n != -1)
            {
                if(findIndex(molIds_, molI->id()) != -1)
                {
                    if(componentControl_)
                    {
                        if(X_)
                        {
                            molI->v().x() += deltaU[n].x();
                        }
                        if(Y_)
                        {
                            molI->v().y() += deltaU[n].y();
                        }
                        if(Z_)
                        {
                            molI->v().z() += deltaU[n].z();
                        }
                    }
                    else
                    {
                        molI->v() += deltaU[n];
                    } 
                }
            }
        }
    }    
   
}

void polyVelocityBins::controlBeforeMove()
{}

void polyVelocityBins::controlBeforeForces()
{}

void polyVelocityBins::controlDuringForces
(
    polyMolecule* molI,
    polyMolecule* molJ
)
{}

void polyVelocityBins::controlAfterForces()
{}


void polyVelocityBins::controlAfterVelocityII()
{}


void polyVelocityBins::calculateProperties()
{}

void polyVelocityBins::output
(
    const fileName& fixedPathName, 
    const fileName& timePath
)
{

}



void polyVelocityBins::updateProperties(const dictionary& newDict)
{
    //- the main controller properties should be updated first
    updateStateControllerProperties(newDict);

    propsDict_ = newDict.subDict(typeName + "Properties");
}




} // End namespace Foam

// ************************************************************************* //
