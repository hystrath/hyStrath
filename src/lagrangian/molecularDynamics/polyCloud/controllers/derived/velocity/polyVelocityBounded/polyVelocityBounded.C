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

#include "polyVelocityBounded.H"
#include "addToRunTimeSelectionTable.H"
#include "IFstream.H"
#include "graph.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(polyVelocityBounded, 0);
addToRunTimeSelectionTable(polyStateController, polyVelocityBounded, dictionary);

// * * * * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * * * //

void polyVelocityBounded::setBoundBox()
{
    vector startPoint = propsDict_.lookup("startPoint");
    vector endPoint = propsDict_.lookup("endPoint");

    bb_.resetBoundedBox(startPoint, endPoint);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
polyVelocityBounded::polyVelocityBounded
(
    Time& t,
    polyMoleculeCloud& molCloud,
    const dictionary& dict
)
:
    polyStateController(t, molCloud, dict),
    propsDict_(dict.subDict(typeName + "Properties")),
    molIds_(),
    bb_(),
//     avMass_(0.0),
    lambda_(readScalar(propsDict_.lookup("lambda"))),
    deltaT_(t.deltaT().value()), 
    accumulatedTime_(0.0),
    startTime_(0.0),
    endTime_(GREAT),
    X_(false),
    Y_(false),
    Z_(false)
{
    writeInTimeDir_ = false;
    writeInCase_ = false;

//     singleValueController() = true;

    setBoundBox();
    
    molIds_.clear();

    selectIds ids
    (
        molCloud_.cP(),
        propsDict_
    );

    molIds_ = ids.molIds();

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
                FatalErrorIn("polyVelocityBounded::polyVelocityBounded()")
                    << "At least one component (X, Y, Z) should be chosen " << nl << "in: "
                    << time_.time().system()/"controllersDict"
                    << exit(FatalError);
            }
        }
    }
    
    if (propsDict_.found("startAtTime"))
    {    
        startTime_ = readScalar(propsDict_.lookup("startAtTime"));
    }
    
    if (propsDict_.found("endAtTime"))
    {
        endTime_ = readScalar(propsDict_.lookup("endAtTime"));
    }    
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

polyVelocityBounded::~polyVelocityBounded()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void polyVelocityBounded::initialConfiguration()
{
// 	readProperties();
}



void polyVelocityBounded::controlBeforeVelocityI()
{
    scalar mass = 0.0;
    scalar mols = 0.0;
    vector mom = vector::zero;
    
    {
        IDLList<polyMolecule>::iterator mol(molCloud_.begin());

        for (mol = molCloud_.begin(); mol != molCloud_.end(); ++mol)
        {
            if(findIndex(molIds_, mol().id()) != -1)
            {
                if(bb_.contains(mol().position()))
                {
                    const scalar& massI = molCloud_.cP().mass(mol().id());

                    mass += massI;
                    mom += mol().v()*massI;
                    mols += 1.0;
                }
            }
        }
    }

    if (Pstream::parRun())
    {
        reduce(mols, sumOp<scalar>());        
        reduce(mass, sumOp<scalar>());
        reduce(mom, sumOp<vector>());
    }
    
    vector velocityMeasured = vector::zero;
    
    if(mass > 0)
    {
        velocityMeasured = mom/mass;
    }

    vector deltaU = (velocity_ - velocityMeasured)*lambda_;

    accumulatedTime_ += deltaT_;
    
    if((accumulatedTime_ >= startTime_) && (accumulatedTime_ <= endTime_))
    {    
        IDLList<polyMolecule>::iterator mol(molCloud_.begin());

        for (mol = molCloud_.begin(); mol != molCloud_.end(); ++mol)
        {
            if(findIndex(molIds_, mol().id()) != -1)
            {
                if(bb_.contains(mol().position()))
                {
                    if(componentControl_)
                    {
                        if(X_)
                        {
                            mol().v().x() += deltaU.x();
                        }
                        if(Y_)
                        {
                            mol().v().y() += deltaU.y();
                        }
                        if(Z_)
                        {
                            mol().v().z() += deltaU.z();
                        }
                    }
                    else
                    {
                        mol().v() += deltaU;
                    }
                }
            }
        }
    
        Info << "polyVelocityBounded - controlling: " << " target velocity = " << velocity_ 
                << ", vel Meas = " << velocityMeasured  
                << ", DeltaU = " << deltaU
                << endl;
    }
}

void polyVelocityBounded::controlBeforeMove()
{}

void polyVelocityBounded::controlBeforeForces()
{}

void polyVelocityBounded::controlDuringForces
(
    polyMolecule* molI,
    polyMolecule* molJ
)
{}

void polyVelocityBounded::controlAfterForces()
{}

void polyVelocityBounded::controlAfterVelocityII()
{}



void polyVelocityBounded::calculateProperties()
{}

void polyVelocityBounded::output
(
    const fileName& fixedPathName,
    const fileName& timePath
)
{}

void polyVelocityBounded::updateProperties(const dictionary& newDict)
{
    //- the main controller properties should be updated first
    updateStateControllerProperties(newDict);

    propsDict_ = newDict.subDict(typeName + "Properties");
}

} // End namespace Foam

// ************************************************************************* //
