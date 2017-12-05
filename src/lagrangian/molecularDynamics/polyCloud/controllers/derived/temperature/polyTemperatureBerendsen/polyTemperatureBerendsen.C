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

#include "polyTemperatureBerendsen.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(polyTemperatureBerendsen, 0);
addToRunTimeSelectionTable(polyStateController, polyTemperatureBerendsen, dictionary);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
polyTemperatureBerendsen::polyTemperatureBerendsen
(
    Time& t,
    polyMoleculeCloud& molCloud,
    const dictionary& dict
)
:
    polyStateController(t, molCloud, dict),
    propsDict_(dict.subDict(typeName + "Properties")),
 	tauT_(readScalar(propsDict_.lookup("tauT"))),
    componentControl_(false),
    X_(false),
    Y_(false),
    Z_(false),
    peculiar_(false)
{
    temperature_ = readScalar(propsDict_.lookup("temperature"));
    
    molIds_.clear();

    selectIds ids
    (
        molCloud_.cP(),
        propsDict_
    );

    molIds_ = ids.molIds();    
    
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
                FatalErrorIn("polyTemperatureBerendsen::polyTemperatureBerendsen()")
                    << "At least one component (X, Y, Z) should be chosen " << nl << "in: "
                    << time_.time().system()/"controllersDict"
                    << exit(FatalError);
            }
        }
    }

    if (propsDict_.found("peculiar"))
    {
        peculiar_ = Switch(propsDict_.lookup("peculiar"));
    }
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

polyTemperatureBerendsen::~polyTemperatureBerendsen()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void polyTemperatureBerendsen::initialConfiguration()
{}

void polyTemperatureBerendsen::controlBeforeVelocityI()
{}

void polyTemperatureBerendsen::controlBeforeMove()
{}

void polyTemperatureBerendsen::controlBeforeForces()
{}

void polyTemperatureBerendsen::controlDuringForces
(
    polyMolecule* molI,
    polyMolecule* molJ
)
{}

void polyTemperatureBerendsen::controlAfterForces()
{}

void polyTemperatureBerendsen::controlAfterVelocityII()
{
	const scalar deltaTMD = time_.deltaT().value(); 
    
    if(control_)
    {
        // - calculate streaming velocity
        scalar mass = 0.0;
        vector mom = vector::zero;

        forAll(controlZone(), c)
        {
            const label& cellI = controlZone()[c];
            const List<polyMolecule*>& molsInCell = molCloud_.cellOccupancy()[cellI];
        
            forAll(molsInCell, m)
            {
                polyMolecule* molI = molsInCell[m];
                
                if(findIndex(molIds_, molI->id()) != -1)
                {
                    const scalar& massI = molCloud_.cP().mass(molI->id());

                    mass += massI;
                    mom += massI*molI->v();
                }
            }
        }

        if (Pstream::parRun())
        {
            reduce(mom, sumOp<vector>());
            reduce(mass, sumOp<scalar>());
        }

        vector velocity(vector::zero);

        if(mass > 0.0)
        {
            velocity = mom/mass;
        }

        scalar kE = 0.0;
        scalar dof = 0.0;
        scalar angularKeSum = 0.0;

        // - kinetic energy (from thermal velocities)
        forAll(controlZone(), c)
        {
            const label& cellI = controlZone()[c];
            const List<polyMolecule*>& molsInCell = molCloud_.cellOccupancy()[cellI];
        
            forAll(molsInCell, m)
            {
                polyMolecule* molI = molsInCell[m];


                if(findIndex(molIds_, molI->id()) != -1)            
                {
                    const scalar& massI = molCloud_.cP().mass(molI->id());

                    kE += 0.5*massI*magSqr(molI->v() - velocity);
                    dof += molCloud_.cP().degreesOfFreedom(molI->id());

                    const diagTensor& molMoI(molCloud_.cP().momentOfInertia(molI->id()));

                    // angular speed 
                    const vector& molOmega(inv(molMoI) & molI->pi());
                    angularKeSum += 0.5*(molOmega & molMoI & molOmega);
                }
            }
        }

        if (Pstream::parRun())
        {
            reduce(kE, sumOp<scalar>());
            reduce(dof, sumOp<scalar>());
            reduce(angularKeSum, sumOp<scalar>());
        }

        scalar tempMeasI = 0.0;

        if(dof > 0.0)
        {
            const scalar& kB = molCloud_.redUnits().kB();

            tempMeasI = (2.0*(kE+angularKeSum))/(kB*dof);

            const reducedUnits& rU = molCloud_.redUnits();

            Info<< "Temp Berendsen, zone : " << regionName() 
                << " T = " << tempMeasI << " (reduced units) "
                << " T = " << tempMeasI*rU.refTemp() << " (SI units) "
                << endl;
        }
        else
        {
            tempMeasI = temperature_;

            Info << "WARNING: dof is zero" << endl;
        }

        scalar chi = sqrt(1.0 + (deltaTMD/tauT_)*((temperature_/tempMeasI) - 1.0) );

        forAll(controlZone(), c)
        {
            const label& cellI = controlZone()[c];
            const List<polyMolecule*>& molsInCell = molCloud_.cellOccupancy()[cellI];

            forAll(molsInCell, m)
            {
                polyMolecule* molI = molsInCell[m];
                
                if(findIndex(molIds_, molI->id()) != -1)            
                {
                    if(!peculiar_)
                    {
                        if(componentControl_)
                        {
                            if(X_)
                            {
                                molI->v().x() *= chi;
                            }
                            if(Y_)
                            {
                                molI->v().y() *= chi;
                            }
                            if(Z_)
                            {
                                molI->v().z() *= chi;
                            }
                        }
                        else
                        {
                            molI->v() *= chi;
                        }
                    }
                    else
                    {
                        molI->v() = velocity + (chi*(molI->v() - velocity));
                    }
                }
            }
        }
    }
}

void polyTemperatureBerendsen::calculateProperties()
{}

void polyTemperatureBerendsen::output
(
    const fileName& fixedPathName, 
    const fileName& timePath
)
{}

void polyTemperatureBerendsen::updateProperties(const dictionary& newDict)
{
    //- the main controller properties should be updated first
    updateStateControllerProperties(newDict);

    propsDict_ = newDict.subDict(typeName + "Properties");

//     readProperties();
}

} // End namespace Foam

// ************************************************************************* //
