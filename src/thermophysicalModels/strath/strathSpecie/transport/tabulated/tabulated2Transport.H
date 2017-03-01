/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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
    Foam::tabulated2Transport

Description
    Transport package using interpolation tables.

    Templated into a given thermodynamics package (needed for thermal
    conductivity).

SourceFiles
    tabulated2TransportI.H
    tabulated2Transport.C

\*---------------------------------------------------------------------------*/

#ifndef tabulated2Transport_H
#define tabulated2Transport_H

#include "interpolation2DTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

template<class Thermo> class tabulated2Transport;

template<class Thermo>
inline tabulated2Transport<Thermo> operator+
(
    const tabulated2Transport<Thermo>&,
    const tabulated2Transport<Thermo>&
);

template<class Thermo>
inline tabulated2Transport<Thermo> operator-
(
    const tabulated2Transport<Thermo>&,
    const tabulated2Transport<Thermo>&
);

template<class Thermo>
inline tabulated2Transport<Thermo> operator*
(
    const scalar,
    const tabulated2Transport<Thermo>&
);

template<class Thermo>
inline tabulated2Transport<Thermo> operator==
(
    const tabulated2Transport<Thermo>&,
    const tabulated2Transport<Thermo>&
);

template<class Thermo>
Ostream& operator<<
(
    Ostream&,
    const tabulated2Transport<Thermo>&
);


/*---------------------------------------------------------------------------*\
                    Class tabulated2Transport Declaration
\*---------------------------------------------------------------------------*/

template<class Thermo>
class tabulated2Transport
:
    public Thermo
{
    // Private data

        //- dynamic viscosity [kg/ms]
        interpolation2DTable<scalar> mu_;

	      // thermal conductivity [W/(Km)]
        interpolation2DTable<scalar> kappa_;


    // Private Member Functions

        //- Construct from components
        inline tabulated2Transport
        (
            const Thermo& t,
            const interpolation2DTable<scalar> mu,
            const interpolation2DTable<scalar> kappa
        );


public:

    // Constructors

        //- Construct as named copy
        inline tabulated2Transport(const word&, const tabulated2Transport&);

        //- Construct from Istream
        //tabulated2Transport(Istream&);

        //- Construct from dictionary
        tabulated2Transport(const dictionary& dict);

        //- Construct and return a clone
        inline autoPtr<tabulated2Transport> clone() const;

        // Selector from Istream
        inline static autoPtr<tabulated2Transport> New(Istream& is);

        // Selector from dictionary
        inline static autoPtr<tabulated2Transport> New(const dictionary& dict);


    // Member functions

        //- Return the instantiated type name
        static word typeName()
        {
            return "tabulated<" + Thermo::typeName() + '>';
        }

        //- Dynamic viscosity [kg/ms]
        inline scalar mu(const scalar p, const scalar T) const;

        //- Thermal conductivity [W/mK]
        inline scalar kappa(const scalar p, const scalar T) const;

        //- Write to Ostream
        void write(Ostream& os) const;


    // Member operators

        inline tabulated2Transport& operator=(const tabulated2Transport&);

        inline void operator+=(const tabulated2Transport&);

        inline void operator-=(const tabulated2Transport&);

        inline void operator*=(const scalar);


    // Friend operators

        friend tabulated2Transport operator+ <Thermo>
        (
            const tabulated2Transport&,
            const tabulated2Transport&
        );

        friend tabulated2Transport operator- <Thermo>
        (
            const tabulated2Transport&,
            const tabulated2Transport&
        );

        friend tabulated2Transport operator* <Thermo>
        (
            const scalar,
            const tabulated2Transport&
        );

        friend tabulated2Transport operator== <Thermo>
        (
            const tabulated2Transport&,
            const tabulated2Transport&
        );


    // Ostream Operator

        friend Ostream& operator<< <Thermo>
        (
            Ostream&,
            const tabulated2Transport&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "tabulated2TransportI.H"

#ifdef NoRepository
#   include "tabulated2Transport.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
