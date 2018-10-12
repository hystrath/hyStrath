/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2005 OpenCFD Ltd.
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
    pdFaceTracker

Description

SourceFiles
    pdFaceTrackerI.H
    pdFaceTracker.C

\*---------------------------------------------------------------------------*/

#ifndef pdFaceTracker_H
#define pdFaceTracker_H

#include "vector.H"
#include "fvMesh.H"
#include "pdParcel.H"
#include "surfaceFields.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class pdFaceTracker Declaration
\*---------------------------------------------------------------------------*/
// class pdCloud;


class pdFaceTracker
{
    // Private data

        const fvMesh& mesh_;

        pdCloud& cloud_;

        List<scalarField> parcelIdFlux_;
        List<scalarField> massIdFlux_;




    // Private Member Functions



public:

    friend class pdParcel;

    // Constructors

        //- Construct from mesh and cloud
        pdFaceTracker
        (
            const polyMesh&,
            pdCloud& cloud
        );

        //- Construct from mesh, cloud and boolean (pdFoam)
        pdFaceTracker
        (
            const polyMesh&,
            pdCloud& cloud,
            const bool& dummy
        );

    // Destructor

        ~pdFaceTracker();

    // Member Functions

        void updateFields
        (
            pdParcel& p
        );

        void clean();

        // Access


            inline const List<scalarField>& parcelIdFlux() const;

            inline const List<scalarField>& massIdFlux() const;


        // Check

        // Edit

        // Write

    // Member Operators

    // Friend Functions

    // Friend Operators


    // IOstream Operators

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "pdFaceTrackerI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
