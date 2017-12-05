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

    Gravity force - gaussian distribution (stress next to a wall)

SourceFiles
    standardGaussianI.H
    standardGaussian.C
    standardGaussianIO.C

\*---------------------------------------------------------------------------*/

#ifndef standardGaussian_H
#define standardGaussian_H

#include "gravityForce.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class standardGaussian Declaration
\*---------------------------------------------------------------------------*/

class standardGaussian
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

        scalar amplitude_;

        scalar constantB_;

        vector forceDirection_;



public:

    //- Runtime type information
    TypeName("standardGaussian");

    // Static data members

    // Constructors

        //- construct from components
        standardGaussian
        (
            Time& time,
            const dictionary& dict
        );

    // Selectors


    // Destructor

        ~standardGaussian();

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
