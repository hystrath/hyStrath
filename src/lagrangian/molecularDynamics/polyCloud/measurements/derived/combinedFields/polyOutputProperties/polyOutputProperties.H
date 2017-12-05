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

Description
    The model which outputs standard properties to the terminal screen 
    in an MD simulation run.

SourceFiles
    polyOutputPropertiesI.H
    polyOutputProperties.C
    polyOutputPropertiesIO.C

\*---------------------------------------------------------------------------*/

#ifndef polyOutputProperties_H
#define polyOutputProperties_H

#include "polyField.H"
#include "vector.H"
#include "polyMesh.H"
#include "volFields.H"
#include "polyMoleculeCloud.H"
// #include "distribution.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class polyOutputProperties Declaration
\*---------------------------------------------------------------------------*/

class polyOutputProperties
:
	public polyField
{

private:

    // Private data

        propertyField fields_;

        vector accumulatedTotalLinearMomentum_;
        scalar accumulatedTotalMass_;
        scalar accumulatedTotalAngularKE_;
        scalar accumulatedTotalLinearKE_;
        scalar accumulatedTotalPE_;
        scalar accumulatedTotalrDotfSum_;
        scalar accumulatedNMols_;
        scalar accumulatedDOFs_;
        scalar averageTemperature_;
        scalar averagePressure_;

        scalar meshVolume_;
        scalar nAvTimeSteps_;
 

    // Private Member Functions

public:

    //- Runtime type information
    TypeName("polyOutputProperties");

    // Constructors

        //- Construct from components
        polyOutputProperties
		(
            Time& t,
            const polyMesh& mesh,
            polyMoleculeCloud& molCloud,
            const dictionary& dict
		);

    // Destructor

 		~polyOutputProperties();

    // Member Functions

        //- create initial field
        void createField();

        //- calculate field
        void calculateField();

        //- write field
        void writeField();

        void measureDuringForceComputation
        (
            polyMolecule* molI,
            polyMolecule* molJ
        );

        void measureDuringForceComputationSite
        (
            polyMolecule* molI,
            polyMolecule* molJ,
            label sI,
            label sJ
        );

        const propertyField& fields() const;
}; 		


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
