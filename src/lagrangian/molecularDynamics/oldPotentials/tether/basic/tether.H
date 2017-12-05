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
    An abstract model for applying gravitational forces

Description

SourceFiles
    tetherI.H
    tether.C
    tetherIO.C

\*---------------------------------------------------------------------------*/

#ifndef tether_H
#define tether_H

#include "IOdictionary.H"
#include "autoPtr.H"
#include "Time.H"
#include "polyMesh.H"
#include "mathematicalConstants.H"
#include "reducedUnits.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class tether Declaration
\*---------------------------------------------------------------------------*/

class tether
{

protected:

        Time& time_;

private:

    // Private Member Functions


public:

    //- Runtime type information
    TypeName("tether");

    // Declare runtime constructor selection table
        declareRunTimeSelectionTable
        (
            autoPtr,
            tether,
            dictionary,
            (
                Time& time,
                const reducedUnits& rU,
                const dictionary& dict
            ),
            (time, rU, dict)
        );

    // Static data members

    // Constructors

        //- Construct from components
        tether
        (
            Time& time,
            const reducedUnits& rU,
            const dictionary& dict
        );


    // Selectors

        static autoPtr<tether> New
        (
            Time& time,
            const reducedUnits& rU,
            const dictionary& dict
        );

    // Destructor

        virtual ~tether();

    // Member Functions


        // Access

            virtual vector force(const vector& rIJ) = 0;
            
            virtual scalar potentialEnergy(const scalar& rIJMag) = 0;
            
//             virtual void updateForce() = 0;

        // Check

        // Edit

        // Write
            virtual void output() = 0;
            virtual void write
            (
                const fileName& fixedPathName,
                const fileName& timePath
            ) = 0;

            virtual void updateProperties
            (
                const dictionary& dict
            ) = 0;


            //- output data

    // Member Operators

    // Friend Functions

    // Friend Operators

    // IOstream Operators

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "tetherI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
