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
    pdPatchBoundary

Description

    Used for applying instantaneous actions on pdParcels colliding with a boundary
    patch.

SourceFiles
    pdPatchBoundaryI.H
    pdPatchBoundary.C
    pdPatchBoundaryIO.C

\*---------------------------------------------------------------------------*/

#ifndef pdPatchBoundary_H
#define pdPatchBoundary_H

#include "IOdictionary.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "vector.H"
#include "volFields.H"
#include "Random.H"
#include "pdParcel.H"
#include "timeData.H"
#include "timeDataMeas.H"
// #include "propertyField.H"
#include "writeTimeData.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class pdPatchBoundary Declaration
\*---------------------------------------------------------------------------*/

class pdPatchBoundary
{

protected:

    // Protected data

        const fvMesh& mesh_;

        pdCloud& cloud_;

        Time& time_;

    // 		dictionary timeDict_;

    //         timeDataMeas time_;

    //         Random rndGen_;

        //- subDictionary containing the properties
        dictionary boundaryDict_;

        //- name of boundary patch
        word patchName_;
        label patchId_;

        //- faces on the boundary
        labelList faces_;
        label nFaces_;
        scalar patchSurfaceArea_;
        scalar totalPatchSurfaceArea_;

        //- cells corresponding to the boundary faces
        labelList cells_;

        //- set all the properties below from model if required

        scalar density_;
        vector velocity_;
        scalar temperature_;

        //- targeted fields
        scalarField densities_;
        vectorField velocities_;
        scalarField temperatures_;

        //- writing out properties
        bool writeInTimeDir_;
        bool writeInCase_;

        bool measurePropertiesAtWall_;

        //- measured fields at the patch
        //scalar preIZ_;
        vector preJP_;
        label  preID_;
        scalar preIE_;
        vector preIMom_;



    // Private Member Functions


public:

    //- Runtime type information
    TypeName("pdPatchBoundary");

    // Declare runtime constructor selection table
        declareRunTimeSelectionTable
        (
            autoPtr,
            pdPatchBoundary,
            dictionary,
            (
                Time& t,
                const polyMesh& mesh,
                pdCloud& molCloud,
                const dictionary& dict
            ),
            (t, mesh, molCloud, dict)
        );

    // Constructors

        //- Construct from components
        pdPatchBoundary
        (
            Time& t,
            const polyMesh& mesh,
            pdCloud& molCloud,
            const dictionary& dict
        );


    // Selectors

        static autoPtr<pdPatchBoundary> New
        (
            Time& t,
            const polyMesh& mesh,
            pdCloud& molCloud,
            const dictionary& dict
        );


    // Destructor

        virtual ~pdPatchBoundary();


    // Member Functions

        void measurePropertiesBeforeControl(pdParcel& p);
        void measurePropertiesAfterControl(pdParcel& p);

        void setBoundaryFields();

//         void writeTimeData
//         (
//             const fileName& pathName,
//             const word& fileName,
//             const List< Pair<scalar> >& data
//         );

        //- scalar field
        void writeTimeData
        (
            const fileName& pathName,
            const word& fileName,
            const scalarField& xAxis
//             const scalarField& yAxis
        );

//         //- vector field
//         void writeTimeData
//         (
//             const fileName& pathName,
//             const word& nameFile,
//             const scalarField& xData,
//             const vectorField& yData
//         );

        //- calculate any required properties
        virtual void calculateProperties() = 0;

        virtual void initialConfiguration() = 0;

        //- control the pdParcel
        virtual void controlParticle
        (
            pdParcel& mol,
            pdParcel::trackingData& td
        ) = 0;

        //- output data
        virtual void output
        (
            const fileName& fixedPathName,
            const fileName& timePath
    ) = 0;

//         virtual void writeField
//         (
//         ) = 0;

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

        //- return the targeted fields
        const scalarField& densityField() const;
        scalarField& densityField();

        const vectorField& velocityField() const;
        vectorField& velocityField();

        const scalarField& temperatureField() const;
        scalarField& temperatureField();

        const bool& writeInTimeDir() const;
        const bool& writeInCase() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
