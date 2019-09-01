/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

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
    Foam::HoffertLienReactionRate

Description
    Hoffert and Lien (1967) reaction rate given by:

        k = A * T^beta * (2 + Ta/T) * exp(-Ta/T)

SourceFiles
    HoffertLienReactionRateI.H

Journal article    
    Quasi-One-Dimensional, Nonequilibrium Gas Dynamics of Partially Ionized 
        Two-Temperature Argon 
    Martin I. Hoffert, and Hwachii Lien
    The Physics of Fluids 10, 1769 (1967)
    doi: 10.1063/1.1762356    

\*---------------------------------------------------------------------------*/

#ifndef HoffertLienReactionRate_H
#define HoffertLienReactionRate_H

#include "scalarField.H"
#include "typeInfo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class HoffertLienReactionRate Declaration
\*---------------------------------------------------------------------------*/

class HoffertLienReactionRate
{
    // Private data

        scalar A_;
        scalar beta_;
        scalar Ta_;


public:

    // Constructors

        //- Construct from components
        inline HoffertLienReactionRate
        (
            const scalar A,
            const scalar beta,
            const scalar Ta
        );

        //- Construct from Istream
        inline HoffertLienReactionRate
        (
            const speciesTable& species,
            Istream& is
        );

        //- Construct from dictionary
        inline HoffertLienReactionRate
        (
            const speciesTable& species,
            const dictionary& dict
        );


    // Member Functions

        //- Return the type name
        static word type()
        {
            return "HoffertLien";
        }

        inline scalar operator()
        (
            const scalar p,
            const scalar T,
            const scalarField& c
        ) const;

        //- Write to stream
        inline void write(Ostream& os) const;


    // Ostream Operator

        inline friend Ostream& operator<<
        (
            Ostream&,
            const HoffertLienReactionRate&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "HoffertLienReactionRateI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
