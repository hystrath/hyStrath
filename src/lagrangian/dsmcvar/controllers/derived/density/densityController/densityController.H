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
    densityController

Description
	This is a density controller, which uses a slightly modified version of the 
    particle insertion protocol called USHER, proposed by Delgado-Buscalioni and Coveney.

SourceFiles
    densityController.C


\*---------------------------------------------------------------------------*/

#ifndef densityController_H
#define densityController_H

#include "dsmcStateController.H"
#include "vector.H"
#include "fvMesh.H"
#include "volFields.H"
#include "dsmcCloud.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class densityController Declaration
\*---------------------------------------------------------------------------*/

class densityController
:
   public dsmcStateController
{

private:
    // Private data

    //- subDictionary containing the properties
    dictionary propsDict_;
    label typeId_;

    //- initial target density, velocity and temperature
    scalar density_; // parcel density
    vector velocity_;
    scalar temperature_;
	scalarField densities_;
	vectorField velocities_;
	scalarField temperatures_;
    labelField nMolsToControl_;
    scalarField parcels_;
    scalarField measuredDensity_;
    scalarField residual_;
    scalarField residualSum_;

    // Private Member Functions
    void deleteParcelFromDSMC(const label& c);
    void insertParcelWithinDSMC(const label& c);
	
	void setProperties();


public:

    //- Runtime type information
    TypeName("densityController");

    // Constructors

    //- Construct from components
    densityController
    (
		Time& t,
        dsmcCloud& dsmcCloud,
        const dictionary& dict
    );

    // Destructor
    ~densityController();

    // Member Functions

    void initialConfiguration();
	
	void controlParcelsBeforeMove();
	void controlParcelsBeforeCollisions();
	void controlParcelsAfterCollisions();

    //- calculate any required properties
    void calculateProperties();

    //- output data
    void output
    (
	const fileName& fixedPathName, 
	const fileName& timePath
    );

    //- update properties from a modified dictionary
    void updateProperties(const dictionary&);

    //- set no of molecules to control
    void nMolsToControl();

        // Access
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
