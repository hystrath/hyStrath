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
    pdGeneralBoundary

Description

SourceFiles
    pdGeneralBoundaryI.H
    pdGeneralBoundary.C
    pdGeneralBoundaryIO.C

\*---------------------------------------------------------------------------*/

#ifndef pdGeneralBoundary_H
#define pdGeneralBoundary_H

#include "IOdictionary.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "vector.H"
#include "volFields.H"
#include "polyMesh.H"
#include "pdParcel.H"
#include "timeFluxData.H"
// #include "timeDataMeas.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class pdGeneralBoundary Declaration
\*---------------------------------------------------------------------------*/

class pdGeneralBoundary
{

protected:

    // Protected data

        const fvMesh& mesh_;

        pdCloud& cloud_;

//         Time& time_;
        dictionary boundaryDict_;

//         dictionary timeDict_;

//         timeFluxData time_;
//         timeDataMeas timeMeas_;

        //- name of boundary patch
        word patchName_;
        label patchId_;

        //- faces on the boundary
        labelList faces_;
        label nFaces_;

        // patch surface area
        scalar patchSurfaceArea_;

        //- cells corresponding to the boundary faces
        labelList cells_;

        //- set all the properties below from model if required

            scalar density_;
            vector velocity_;
            scalar temperature_;
            /**************************************************************/
            //scalar charge_;
            /**************************************************************/

            //- targeted fields
            scalarField densities_;
            vectorField velocities_;
            scalarField temperatures_;

            /**************************************************************/
            //scalarField charge_;
            /**************************************************************/

            //- writing out properties
            bool writeInTimeDir_;
            bool writeInCase_;

//         Random rndGen_;

    // Private Member Functions

public:

    //- Runtime type information
    TypeName("pdGeneralBoundary");

    // Declare runtime constructor selection table
        declareRunTimeSelectionTable
        (
            autoPtr,
            pdGeneralBoundary,
            dictionary,
            (
                Time& t,
                const polyMesh& mesh,
                pdCloud& cloud,
                const dictionary& dict
            ),
            (t, mesh, cloud, dict)
        );

    // Constructors

        //- Construct from components
        pdGeneralBoundary
        (
            Time& t,
            const polyMesh& mesh,
            pdCloud& cloud,
            const dictionary& dict
        );


    // Selectors

        static autoPtr<pdGeneralBoundary> New
        (
            Time& t,
            const polyMesh& mesh,
            pdCloud& cloud,
            const dictionary& dict
        );


    // Destructor

        virtual ~pdGeneralBoundary();


    // Member Functions

        void updateTime();


        //- create an initial configuration
        virtual void initialConfiguration() = 0;

        //- control
 		virtual void controlParcelsBeforeMove() = 0;
        virtual void controlParcelsBeforeCollisions() = 0;
        virtual void controlParcelsAfterCollisions() = 0;

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

        void writeTimeData
        (
            const fileName& pathName,
            const word& nameFile,
            const scalarField& xData,
            const scalarField& yData
        );

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
            const scalar& density() const;
            scalar& density();

            const vector& velocity() const;
            vector& velocity();

            const scalar& temperature() const;
            scalar& temperature();

            /**************************************************************/
            /*
            const scalar& charge() const;
            scalar& charge();
            */
            /**************************************************************/

            //- return the targeted fields
            const scalarField& densityField() const;
            scalarField& densityField();

            const vectorField& velocityField() const;
            vectorField& velocityField();

            const scalarField& temperatureField() const;
            scalarField& temperatureField();

            /**************************************************************/
            /*
            const scalarField& chargeField() const;
            scalarField& chargeField();
            */
            /**************************************************************/

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
