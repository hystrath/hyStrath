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
    An abstract model for moving/controlling MD walls

Description

SourceFiles
    wallMotionI.H
    wallMotion.C
    wallMotionIO.C

\*---------------------------------------------------------------------------*/

#ifndef wallMotion_H
#define wallMotion_H

#include "IOdictionary.H"
#include "autoPtr.H"
#include "Time.H"
#include "polyMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class wallMotion Declaration
\*---------------------------------------------------------------------------*/

class wallMotion
{

protected:

        Time& time_;

private:

    // Private Member Functions


public:

    //- Runtime type information
    TypeName("wallMotion");

    // Declare runtime constructor selection table
        declareRunTimeSelectionTable
        (
            autoPtr,
            wallMotion,
            dictionary,
            (
                Time& time,
                const dictionary& dict
            ),
            (time, dict)
        );

    // Static data members

    // Constructors

        //- Construct from components
        wallMotion
        (
            Time& time,
            const dictionary& dict
        );


    // Selectors

        static autoPtr<wallMotion> New
        (
            Time& time,
            const dictionary& dict
        );

    // Destructor

        virtual ~wallMotion();

    // Member Functions


        // Access

            virtual const vector& velocity() const = 0;
            virtual void updateVelocity() = 0;
        // Check

        // Edit

        // Write

            //- output data

    // Member Operators

    // Friend Functions

    // Friend Operators

    // IOstream Operators

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "wallMotionI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
