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

#include "polyVelocity2DBins.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(polyVelocity2DBins, 0);

addToRunTimeSelectionTable(polyStateController, polyVelocity2DBins, dictionary);



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
polyVelocity2DBins::polyVelocity2DBins
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
    molIds_(),
    deltaT_(t.deltaT().value())
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
                FatalErrorIn("polyVelocity2DBins::polyVelocity2DBins()")
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

    binModel_ =  autoPtr<twoDimBinModel>
    (
        twoDimBinModel::New(mesh_, propsDict_)
    );

}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

polyVelocity2DBins::~polyVelocity2DBins()
{}



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void polyVelocity2DBins::initialConfiguration()
{
}

void polyVelocity2DBins::controlBeforeVelocityI()
{
/*    const List< DynamicList<polyMolecule*> >& cellOccupancy
        = molCloud_.cellOccupancy();

    const labelList& cells = mesh_.cellZones()[regionId_];
    
    const List<label>& noOfBins = binModel_->nBins();
    
    label nBinsX = noOfBins[0];
    label nBinsY = noOfBins[1];
    
    List<scalarField> mass(nBinsX);
    List<vectorField> momentum(nBinsX);
    
    forAll(mass, i)
    {
        mass[i].setSize(nBinsY, 0.0);
        momentum[i].setSize(nBinsY, vector::zero);
    }

    forAll(cells, c)
    {
        const label& cell = cells[c];
        const List<polyMolecule*>& molsInCell = cellOccupancy[cell];

        forAll(molsInCell, mIC)
        {
            polyMolecule* molI = molsInCell[mIC];

            const vector& rI = molI->position();

            List<label> n = binModel_->isPointWithinBin(rI, cell);
            
            label nX = n[0];
            label nY = n[1];
            
            if( (nX + nY) >= 0)
            {
                if(findIndex(molIds_, molI->id()) != -1)
                {
                    const polyMolecule::constantProperties& constProp 
                                        = molCloud_.constProps(molI->id());

                    mass[nX][nY] += constProp.mass();
                    momentum[nX][nY] += constProp.mass()*molI->v();
                }
            }
        }
    }

    if(Pstream::parRun())
    {
        forAll(mass, i)
        {
            forAll(mass[i], j)
            {
                reduce(momentum[i][j], sumOp<vector>());
                reduce(mass[i][j], sumOp<scalar>());
            }
        }
    }

    List<vectorField> velMeasured(nBinsX);
    List<vectorField> deltaU(nBinsX);
    
    forAll(velMeasured, i)
    {
        velMeasured[i].setSize(nBinsY, vector::zero);
        deltaU[i].setSize(nBinsY, vector::zero);
        
        forAll(velMeasured[i], j)
        {
            if(mass[i][j] > 0.0)
            {
                velMeasured[i][j] = momentum[i][j]/mass[i][j];
                deltaU[i][j] = (velocity_ - velMeasured[i][j])*lambda_;    
                Info << "vel = " << velMeasured[i][j] << ", deltaU = " << deltaU[i][j] << endl;
            }
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

            List<label> n = binModel_->isPointWithinBin(rI, cell);
            
            label nX = n[0];
            label nY = n[1];
            
            if( (nX + nY) >= 0)
            {
                if(findIndex(molIds_, molI->id()) != -1)
                {
                    if(componentControl_)
                    {
                        if(X_)
                        {
                            molI->v().x() += deltaU[nX][nY].x();
                        }
                        if(Y_)
                        {
                            molI->v().y() += deltaU[nX][nY].y();
                        }
                        if(Z_)
                        {
                            molI->v().z() += deltaU[nX][nY].z();
                        }
                    }
                    else
                    {
                        molI->v() += deltaU[nX][nY];
                    } 
                }
            }
        }
    }*/    
   
}

void polyVelocity2DBins::controlBeforeMove()
{}

void polyVelocity2DBins::controlBeforeForces()
{}

void polyVelocity2DBins::controlDuringForces
(
    polyMolecule* molI,
    polyMolecule* molJ
)
{}

void polyVelocity2DBins::controlAfterForces()
{
    const List< DynamicList<polyMolecule*> >& cellOccupancy
        = molCloud_.cellOccupancy();

    const labelList& cells = mesh_.cellZones()[regionId_];
    
    const List<label>& noOfBins = binModel_->nBins();
    
    label nBinsX = noOfBins[0];
    label nBinsY = noOfBins[1];

    List<scalarField> mols(nBinsX);    
    List<scalarField> mass(nBinsX);
    List<vectorField> momentum(nBinsX);
    
    forAll(mass, i)
    {
        mols[i].setSize(nBinsY, 0.0);
        mass[i].setSize(nBinsY, 0.0);
        momentum[i].setSize(nBinsY, vector::zero);
    }

    forAll(cells, c)
    {
        const label& cell = cells[c];
        const List<polyMolecule*>& molsInCell = cellOccupancy[cell];

        forAll(molsInCell, mIC)
        {
            polyMolecule* molI = molsInCell[mIC];

            const vector& rI = molI->position();

            List<label> n = binModel_->isPointWithinBin(rI, cell);
            
            label nX = n[0];
            label nY = n[1];
            
            if( (nX + nY) >= 0)
            {
                if(findIndex(molIds_, molI->id()) != -1)
                {
                    const scalar& massI = molCloud_.cP().mass(molI->id());
                                        
                    mols[nX][nY] += 1.0;
                    mass[nX][nY] += massI;
                    momentum[nX][nY] += massI*molI->v();
                }
            }
        }
    }

    if(Pstream::parRun())
    {
        forAll(mass, i)
        {
            forAll(mass[i], j)
            {
                reduce(mols[i][j], sumOp<scalar>());
                reduce(momentum[i][j], sumOp<vector>());
                reduce(mass[i][j], sumOp<scalar>());
            }
        }
    }

    List<vectorField> velMeasured(nBinsX);
    List<vectorField> deltaU(nBinsX);
    
    forAll(velMeasured, i)
    {
        velMeasured[i].setSize(nBinsY, vector::zero);
        deltaU[i].setSize(nBinsY, vector::zero);
        
        forAll(velMeasured[i], j)
        {
            if(mass[i][j] > 0.0)
            {
                velMeasured[i][j] = momentum[i][j]/mass[i][j];
                deltaU[i][j] = (velocity_ - velMeasured[i][j])*lambda_;    
                
                // temp output
                if( (i == 0) && (j == 0) )
                {
                    Info<< "nMols = " <<  mols[i][j]
                        << ", vel = " << velMeasured[i][j]
                        << ", deltaU = " << deltaU[i][j] << endl;
                }
            }
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

            List<label> n = binModel_->isPointWithinBin(rI, cell);
            
            label nX = n[0];
            label nY = n[1];
            
            if( (nX + nY) >= 0)
            {
                if(findIndex(molIds_, molI->id()) != -1)
                {
//                     const polyMolecule::constantProperties& constProp 
//                                         = molCloud_.constProps(molI->id());

//                     const scalar& massI = constProp.mass();
                    
                    if(componentControl_)
                    {
                        if(X_)
                        {
                            molI->a().x() += deltaU[nX][nY].x()/deltaT_;
                        }
                        if(Y_)
                        {
                            molI->a().y() += deltaU[nX][nY].y()/deltaT_;
                        }
                        if(Z_)
                        {
                            molI->a().z() += deltaU[nX][nY].z()/deltaT_;
                        }
                    }
                    else
                    {
                        molI->a() += deltaU[nX][nY]/deltaT_;
                    } 
                }
            }
        }
    }    
    
}

void polyVelocity2DBins::controlAfterVelocityII()
{}


void polyVelocity2DBins::calculateProperties()
{}

void polyVelocity2DBins::output
(
    const fileName& fixedPathName, 
    const fileName& timePath
)
{}

void polyVelocity2DBins::updateProperties(const dictionary& newDict)
{
    //- the main controller properties should be updated first
    updateStateControllerProperties(newDict);

    propsDict_ = newDict.subDict(typeName + "Properties");
}




} // End namespace Foam

// ************************************************************************* //
