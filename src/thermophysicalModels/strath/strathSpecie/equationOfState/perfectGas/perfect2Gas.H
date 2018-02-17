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
    Foam::perfect2Gas

Description
    Perfect gas equation of state.

SourceFiles
    perfect2GasI.H
    perfect2Gas.C

\*---------------------------------------------------------------------------*/

#ifndef perfect2Gas_H
#define perfect2Gas_H

#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

template<class Specie> class perfect2Gas;

template<class Specie>
inline perfect2Gas<Specie> operator+
(
    const perfect2Gas<Specie>&,
    const perfect2Gas<Specie>&
);

template<class Specie>
inline perfect2Gas<Specie> operator-
(
    const perfect2Gas<Specie>&,
    const perfect2Gas<Specie>&
);

template<class Specie>
inline perfect2Gas<Specie> operator*
(
    const scalar,
    const perfect2Gas<Specie>&
);

template<class Specie>
inline perfect2Gas<Specie> operator==
(
    const perfect2Gas<Specie>&,
    const perfect2Gas<Specie>&
);

template<class Specie>
Ostream& operator<<
(
    Ostream&,
    const perfect2Gas<Specie>&
);


/*---------------------------------------------------------------------------*\
                           Class perfect2Gas Declaration
\*---------------------------------------------------------------------------*/

template<class Specie>
class perfect2Gas
:
    public Specie
{

public:

    // Constructors

        //- Construct from components
        inline perfect2Gas(const Specie& sp);

        //- Construct from Istream
        perfect2Gas(Istream&);

        //- Construct from dictionary
        perfect2Gas(const dictionary& dict);

        //- Construct as named copy
        inline perfect2Gas(const word& name, const perfect2Gas&);

        //- Construct and return a clone
        inline autoPtr<perfect2Gas> clone() const;

        // Selector from Istream
        inline static autoPtr<perfect2Gas> New(Istream& is);

        // Selector from dictionary
        inline static autoPtr<perfect2Gas> New(const dictionary& dict);


    // Member functions

        //- Return the instantiated type name
        static word typeName()
        {
            return "perfect2Gas<" + word(Specie::typeName_()) + '>';
        }


        // Fundamental properties

            //- Is the equation of state is incompressible i.e. rho != f(p)
            static const bool incompressible = false;

            //- Is the equation of state is isochoric i.e. rho = const
            static const bool isochoric = false;

            //- Return density [kg/m^3]
            inline scalar rho(scalar p, scalar T) const;

            //- Return compressibility rho/p [s^2/m^2]
            inline scalar psi(scalar p, scalar T) const;

            //- Return compression factor []
            inline scalar Z(scalar p, scalar T) const;

            //- Return (cp - cv) [J/(kmol K)]
            inline scalar cpMcv(scalar p, scalar T) const;


        // IO

            //- Write to Ostream
            void write(Ostream& os) const;


    // Member operators

        inline void operator+=(const perfect2Gas&);
        inline void operator-=(const perfect2Gas&);

        inline void operator*=(const scalar);


    // Friend operators

        friend perfect2Gas operator+ <Specie>
        (
            const perfect2Gas&,
            const perfect2Gas&
        );

        friend perfect2Gas operator- <Specie>
        (
            const perfect2Gas&,
            const perfect2Gas&
        );

        friend perfect2Gas operator* <Specie>
        (
            const scalar s,
            const perfect2Gas&
        );

        friend perfect2Gas operator== <Specie>
        (
            const perfect2Gas&,
            const perfect2Gas&
        );


    // Ostream Operator

        friend Ostream& operator<< <Specie>
        (
            Ostream&,
            const perfect2Gas&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "perfect2GasI.H"

#ifdef NoRepository
#   include "perfect2Gas.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
