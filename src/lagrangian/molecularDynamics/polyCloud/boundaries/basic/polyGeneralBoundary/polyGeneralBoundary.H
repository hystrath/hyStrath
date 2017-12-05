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
    polyGeneralBoundary

Description

SourceFiles
    polyGeneralBoundaryI.H
    polyGeneralBoundary.C
    polyGeneralBoundaryIO.C

\*---------------------------------------------------------------------------*/

#ifndef polyGeneralBoundary_H
#define polyGeneralBoundary_H

#include "IOdictionary.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "vector.H"
#include "volFields.H"
#include "polyMesh.H"
#include "polyMolecule.H"
#include "writeTimeData.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class polyGeneralBoundary Declaration
\*---------------------------------------------------------------------------*/

class polyGeneralBoundary
{

protected:

    // Protected data

        const fvMesh& mesh_;

        polyMoleculeCloud& molCloud_;

        dictionary boundaryDict_;

        Time& time_; 

        //- name of boundary patch
        word patchName_;
        label patchId_;

        //- faces on the boundary
        labelList faces_;
        label nFaces_;

        //- cells corresponding to the boundary faces
        labelList cells_;

        //- set all the properties below from model if required

		//- targeted fields
		scalarField densities_;
		vectorField velocities_;
		scalarField temperatures_;

		//- writing out properties
		bool writeInTimeDir_;
		bool writeInCase_;

    // Private Member Functions

public:

    //- Runtime type information
    TypeName("polyGeneralBoundary");

    // Declare runtime constructor selection table
        declareRunTimeSelectionTable
        (
            autoPtr,
            polyGeneralBoundary,
            dictionary,
            (
                Time& t,
                const polyMesh& mesh,
                polyMoleculeCloud& molCloud,
                const dictionary& dict
            ),
            (t, mesh, molCloud, dict)
        );
   
    // Constructors

        //- Construct from components
        polyGeneralBoundary
        (
            Time& t,
            const polyMesh& mesh,
            polyMoleculeCloud& molCloud,
            const dictionary& dict
        );


    // Selectors

        static autoPtr<polyGeneralBoundary> New
        (
            Time& t,
            const polyMesh& mesh,
            polyMoleculeCloud& molCloud,
            const dictionary& dict
        );


    // Destructor

        virtual ~polyGeneralBoundary();


    // Member Functions


        //- create an initial configuration
        virtual void initialConfiguration() = 0;

        //- control polyMolecules
 		virtual void controlMols() = 0;

        //- output data
        virtual void output
        (
            const fileName& fixedPathName,
            const fileName& timePath
        ) = 0;

        //- calculate any required properties
        virtual void calculateProperties() = 0;

//         virtual vector computeForce(const vector&) = 0;
//         virtual scalar computeEnergy(const vector&) = 0;
//         virtual scalar computeVirial(const vector&) = 0;


//         virtual void resetSphere() = 0;


        //- update properties from a modified dictionary
        virtual void updateProperties(const dictionary&) = 0;

        void updateBoundaryProperties(const dictionary&);


        // Access

            const labelList& controlPatch() const;

            //- return the adjacent cells
            const labelList& controlZone() const;

            //- return the name of the control patch
            const word& patchName() const;
            const label& patchId() const;


            //- return the targeted fields 
            const scalarField& densityField() const;
            scalarField& densityField();

            const vectorField& velocityField() const;
            vectorField& velocityField();

            const scalarField& temperatureField() const;
            scalarField& temperatureField();

//             scalar avReqDensity();
//             vector avReqVelocity();
//             scalar avReqTemperature();

            const bool& writeInTimeDir() const;
            const bool& writeInCase() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
