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
    dsmcFaceTracker

Description

SourceFiles
    dsmcFaceTrackerI.H
    dsmcFaceTracker.C

\*---------------------------------------------------------------------------*/

#ifndef dsmcFaceTracker_H
#define dsmcFaceTracker_H

#include "vector.H"
#include "fvMesh.H"
#include "dsmcParcel.H"
#include "surfaceFields.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class dsmcFaceTracker Declaration
\*---------------------------------------------------------------------------*/
// class dsmcCloud;


class dsmcFaceTracker
{
    // Private data

        const fvMesh& mesh_;

        dsmcCloud& cloud_;

        List<scalarField> parcelIdFlux_;
        List<scalarField> massIdFlux_;

    // Private Member Functions

public:

    friend class dsmcParcel;

    // Constructors

        //- Construct from mesh and cloud
        dsmcFaceTracker
        (
            const polyMesh&,
            dsmcCloud& cloud
        );

        //- Construct from mesh, cloud and boolean (dsmcFoam)
        dsmcFaceTracker
        (
            const polyMesh&,
            dsmcCloud& cloud,
            const bool& dummy
        );

    // Destructor

        ~dsmcFaceTracker();

    // Member Functions

        void trackParcelFaceTransition
        (
            const dsmcParcel& p
        );

        void trackFaceTransition
        (
            const label& typeId,
            const vector& U,
            const scalar& RWF,
            const label& crossedFaceI
        );

        void clean();

        void reset();

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

#include "dsmcFaceTrackerI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
