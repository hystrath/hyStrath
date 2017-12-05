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

Description

    Gravity force - read in force distribution from file

SourceFiles
    oneDimForceI.H
    oneDimForce.C
    oneDimForceIO.C

\*---------------------------------------------------------------------------*/

#ifndef oneDimForce_H
#define oneDimForce_H

#include "gravityForce.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class oneDimForce Declaration
\*---------------------------------------------------------------------------*/

class oneDimForce
:
    public gravityForce
{

private:

    // Private data

        dictionary propsDict_;

//         vector force_;
        vector startPoint_;
        vector endPoint_;
        vector normalVector_;
        scalar length_;
        vector forceDirection_;

        label nBins_;
        scalar binWidth_;

//         scalarField r_;
        scalarField forces_;


public:

    //- Runtime type information
    TypeName("oneDimForce");

    // Static data members

    // Constructors

        //- construct from components
        oneDimForce
        (
            Time& time,
            const dictionary& dict
        );

    // Selectors


    // Destructor

        ~oneDimForce();

    // Member Functions

        // Access

        // Check

        // Edit

            void updateForce();

            vector force(const vector& position);
            
            vector force(const scalar& time);
            
        // Write


            void write
            (
                const fileName& fixedPathName,
                const fileName& timePath
            );

            scalar returnForce
            (
                const scalar& y
            );

            void updateProperties
            (
                const dictionary& dict
            );

    // Member Operators


    // Friend Functions

    // Friend Operators

    // IOstream Operators
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
