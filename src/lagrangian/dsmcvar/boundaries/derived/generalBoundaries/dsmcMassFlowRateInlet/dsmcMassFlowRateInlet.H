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

Class
    dsmcMassFlowRateInlet

Description

    User defines a fixed mass flow rate and a temperature. The density and velocity
    are calculated as the simulation proceeds. This is the non-parallel communication
    version proposed by Erin Farbar in the "Subsonic Flow Boundary Conditions for the Direct
    Simulation Monte Carlo Method" paper.

    Should be used in conjuction with a deletion patch.

SourceFiles
    dsmcMassFlowRateInlet.C


\*---------------------------------------------------------------------------*/

#ifndef dsmcMassFlowRateInlet_H
#define dsmcMassFlowRateInlet_H

#include "dsmcGeneralBoundary.H"
#include "vector.H"
#include "polyMesh.H"
#include "volFields.H"
#include "dsmcCloud.H"
#include "distribution.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class dsmcMassFlowRateInlet Declaration
\*---------------------------------------------------------------------------*/

class dsmcMassFlowRateInlet
:
	public dsmcGeneralBoundary
{

private:

    // Private data
    dictionary propsDict_;
    labelField typeIds_;
    scalarField moleFractions_;
    vectorField inletVelocity_;
    vectorField previousInletVelocity_;
    vectorField momentum_;
    List <scalarField> n_; // inlet number densites
    List <scalarField> accumulatedParcelsToInsert_;
    scalarField mass_;
    scalar massFlowRate_;
    scalar inletTemperature_;
//     scalar theta_;
    vector initialVelocity_;    

    void setProperties();

public:

    //- Runtime type information
    TypeName("dsmcMassFlowRateInlet");

    // Constructors

        //- Construct from components
    dsmcMassFlowRateInlet
    (
        Time& t,
        const polyMesh& mesh,
        dsmcCloud& cloud,
        const dictionary& dict
    );

    // Destructor

    ~dsmcMassFlowRateInlet();


    // Member Functions

    void initialConfiguration();

    //- calculate any required properties
    void calculateProperties();

    //- control the mols
    void controlParcelsBeforeMove();
    void controlParcelsBeforeCollisions();
    void controlParcelsAfterCollisions();

    //- output data
    void output
    (
        const fileName& fixedPathName,
        const fileName& timePath
    );

    //- update properties from a modified dictionary
    void updateProperties(const dictionary&);
    
    //reset fields if mesh has changed
    void setNewBoundaryFields();

        // Access

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
