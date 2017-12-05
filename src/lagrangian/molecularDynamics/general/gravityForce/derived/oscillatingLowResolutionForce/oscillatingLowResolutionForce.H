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
    Wall is moved in an oscillatory manner.


Description

SourceFiles
    oscillatingLowResolutionForceI.H
    oscillatingLowResolutionForce.C
    oscillatingLowResolutionForceIO.C

\*---------------------------------------------------------------------------*/

#ifndef oscillatingLowResolutionForce_H
#define oscillatingLowResolutionForce_H


#include "gravityForce.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class oscillatingLowResolutionForce Declaration
\*---------------------------------------------------------------------------*/

class oscillatingLowResolutionForce
:
    public gravityForce
{

private:

    // Private data

        dictionary propsDict_;

        vector force_;

        vector unitVector_; 

        scalar omega_;

        scalar amplitude_;

        scalar m_;

//         scalar offsetTime_;

        scalar currentTime_;

        scalar deltaTMD_;

public:

    //- Runtime type information
    TypeName("oscillatingLowResolutionForce");

    // Static data members

    // Constructors

        //- construct from components
        oscillatingLowResolutionForce
        (
            Time& time,
            const dictionary& dict
        );

    // Selectors


    // Destructor

        ~oscillatingLowResolutionForce();

    // Member Functions

        //- output data

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
