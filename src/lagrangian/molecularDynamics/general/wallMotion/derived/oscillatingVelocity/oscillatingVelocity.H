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
    oscillatingVelocityI.H
    oscillatingVelocity.C
    oscillatingVelocityIO.C

\*---------------------------------------------------------------------------*/

#ifndef oscillatingVelocity_H
#define oscillatingVelocity_H


#include "wallMotion.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class oscillatingVelocity Declaration
\*---------------------------------------------------------------------------*/

class oscillatingVelocity
:
    public wallMotion
{

private:

    // Private data

        dictionary propsDict_;

        vector velocity_;
        
        vector unitVector_;
        
        scalar initialVelocity_;
        
        scalar amplitude_;

        scalar period_;

        scalar offsetTime_;

        scalar currentTime_;

//         scalar deltaTMD_;
        
        scalar deltaT_;
        
public:

    //- Runtime type information
    TypeName("oscillatingVelocity");

    // Static data members

    // Constructors

        //- construct from components
        oscillatingVelocity
        (
            Time& time,
            const dictionary& dict
        );

    // Selectors


    // Destructor

        ~oscillatingVelocity();

    // Member Functions

        //- output data

        // Access

        // Check

        // Edit

            void updateVelocity();

            const vector& velocity() const;

        // Write





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
