/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016-2020 hyStrath
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of hyStrath, a derivative work of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    cellMeasurements

Description

SourceFiles
    cellMeasurementsI.H
    cellMeasurements.C

\*---------------------------------------------------------------------------*/

#ifndef cellMeasurements_H
#define cellMeasurements_H

#include "vector.H"
#include "fvMesh.H"
#include "dsmcParcel.H"
#include "surfaceFields.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class cellMeasurements Declaration
\*---------------------------------------------------------------------------*/
// class dsmcCloud;


class cellMeasurements
{
    // Private data

        const fvMesh& mesh_;

        dsmcCloud& cloud_;

        scalarField collisionSeparation_;
        scalarField nColls_;

    // Private Member Functions



public:

    friend class dsmcParcel;

    // Constructors

        //- Construct from mesh and cloud
        cellMeasurements
        (
            const polyMesh&,
            dsmcCloud& cloud
        );

        //- Construct from mesh, cloud and boolean (dsmcFoam)
        cellMeasurements
        (
            const polyMesh&,
            dsmcCloud& cloud,
            const bool& dummy
        );

    // Destructor

        ~cellMeasurements();

    // Member Functions

        void updateFields
        (
            dsmcParcel& p
        );

        void clean();

        void reset();

        // Access
        inline const scalarField& collisionSeparation() const;
        inline const scalarField& nColls() const;

        // Check

        // Edit

        inline scalarField& collisionSeparation();
        inline scalarField& nColls();

        // Write

    // Member Operators

    // Friend Functions

    // Friend Operators


    // IOstream Operators

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "cellMeasurementsI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
