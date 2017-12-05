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

#include "polyVelocityZoneControl.H"
#include "addToRunTimeSelectionTable.H"
#include "IFstream.H"
#include "graph.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(polyVelocityZoneControl, 0);
addToRunTimeSelectionTable(polyStateController, polyVelocityZoneControl, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
polyVelocityZoneControl::polyVelocityZoneControl
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
    Z_(false)
{
    writeInTimeDir_ = false;
    writeInCase_ = false;

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
                FatalErrorIn("polyVelocityZoneControl::polyVelocityZoneControl()")
                    << "At least one component (X, Y, Z) should be chosen " << nl << "in: "
                    << time_.time().system()/"controllersDict"
                    << exit(FatalError);
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

polyVelocityZoneControl::~polyVelocityZoneControl()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void polyVelocityZoneControl::initialConfiguration()
{}

void polyVelocityZoneControl::controlBeforeVelocityI()
{}

void polyVelocityZoneControl::controlBeforeMove()
{}

void polyVelocityZoneControl::controlBeforeForces()
{}

void polyVelocityZoneControl::controlDuringForces
(
    polyMolecule* molI,
    polyMolecule* molJ
)
{}

void polyVelocityZoneControl::controlAfterForces()
{
    // - if control switch is on
	if(control_)
    {
        scalar mass = 0.0;
        vector mom = vector::zero;
        vector velocityMeasured = vector::zero;

        forAll(controlZone(), c)
        {
            const label& cellI = controlZone()[c];
            const List<polyMolecule*>& molsInCell = molCloud_.cellOccupancy()[cellI];

            forAll(molsInCell, m)
            {
                polyMolecule* molI = molsInCell[m];

//                 const polyMolecule::constantProperties& constProp = molCloud_.constProps(molI->id());
//                 const scalar& massI = constProp.mass();
                const scalar& massI = molCloud_.cP().mass(molI->id());
                mass += massI;
                mom += molI->v()*massI;
            }
        }

        if (Pstream::parRun())
        {
            reduce(mass, sumOp<scalar>());
            reduce(mom, sumOp<vector>());
        }

        if(mass > 0)
        {
            velocityMeasured = mom/mass;
        }

        vector deltaU = (velocity_ - velocityMeasured)*lambda_;

        //- rescale the velocities

        const scalar deltaTMD = time_.deltaT().value();

        Info << "polyVelocityZoneControl: " << " target velocity: " << velocity_ 
             << ", vel Meas: " << velocityMeasured << ", delta U: " << deltaU 
             << ", acc: " <<  deltaU/deltaTMD 
             << endl;

        forAll(controlZone(), c)
        {
            const label& cellI = controlZone()[c];
    
            const List<polyMolecule*>& molsInCell = molCloud_.cellOccupancy()[cellI];
    
            forAll(molsInCell, m)
            {
                polyMolecule* molI = molsInCell[m];

                if(componentControl_)
                {
                    if(X_)
                    {
                        molI->a().x() += deltaU.x()/deltaTMD;
                    }
                    if(Y_)
                    {
                        molI->a().y() += deltaU.y()/deltaTMD;
                    }
                    if(Z_)
                    {
                        molI->a().z() += deltaU.z()/deltaTMD;
                    }
                }
                else
                {
                    molI->a() += deltaU/deltaTMD;
                }
            }
        }
    }
}


void polyVelocityZoneControl::controlAfterVelocityII()
{}

void polyVelocityZoneControl::calculateProperties()
{}

void polyVelocityZoneControl::output
(
    const fileName& fixedPathName,
    const fileName& timePath
)
{}

void polyVelocityZoneControl::updateProperties(const dictionary& newDict)
{
    //- the main controller properties should be updated first
    updateStateControllerProperties(newDict);

    propsDict_ = newDict.subDict(typeName + "Properties");

//     lambda_ = readScalar(propsDict_.lookup("lambda"));

//     velocity_ = propsDict_.lookup("velocity");
}

} // End namespace Foam

// ************************************************************************* //
