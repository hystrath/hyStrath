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
    Foam::constantTransport

Description
    Constant properties Transport package for the 2T model.
    Templated into a given thermodynamics package (needed for thermal
    conductivity).

SourceFiles
    constantTransportI.H
    constantTransport.C

\*---------------------------------------------------------------------------*/

#ifndef constantTransport_H
#define constantTransport_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

template<class Thermo> class constantTransport;

template<class Thermo>
inline constantTransport<Thermo> operator+
(
    const constantTransport<Thermo>&,
    const constantTransport<Thermo>&
);

template<class Thermo>
inline constantTransport<Thermo> operator-
(
    const constantTransport<Thermo>&,
    const constantTransport<Thermo>&
);

template<class Thermo>
inline constantTransport<Thermo> operator*
(
    const scalar,
    const constantTransport<Thermo>&
);

template<class Thermo>
inline constantTransport<Thermo> operator==
(
    const constantTransport<Thermo>&,
    const constantTransport<Thermo>&
);

template<class Thermo>
Ostream& operator<<
(
    Ostream&,
    const constantTransport<Thermo>&
);


/*---------------------------------------------------------------------------*\
                           Class constantTransport Declaration
\*---------------------------------------------------------------------------*/

template<class Thermo>
class constantTransport
:
    public Thermo
{
    // Private data

        //- Constant dynamic viscosity [Pa.s]
        scalar mu_;

        //- Constant thermal conductivity, trans-rotational mode [W/mK]
        scalar kappa_;
        
        //- Constant thermal conductivity, vibro-electronic mode [W/mK]
        scalar kappave_;


    // Private Member Functions

        //- Construct from components
        inline constantTransport
        (
            const Thermo& t,
            const scalar mu,
            const scalar kappa,
            const scalar kappave
        );


public:

    // Constructors

        //- Construct as named copy
        inline constantTransport(const word&, const constantTransport&);

        //- Construct from Istream
        constantTransport(Istream&);

        //- Construct from dictionary
        constantTransport(const dictionary& dict);

        //- Construct and return a clone
        inline autoPtr<constantTransport> clone() const;

        // Selector from Istream
        inline static autoPtr<constantTransport> New(Istream& is);

        // Selector from dictionary
        inline static autoPtr<constantTransport> New(const dictionary& dict);


    // Member functions

        //- Return the instantiated type name
        static word typeName()
        {
            return "constant<" + Thermo::typeName() + '>';
        }

        //- Dynamic viscosity [kg/ms]
        inline scalar mu(const scalar p, const scalar T) const;

        //- Thermal conductivity, translational mode [W/mK]
        inline scalar kappatrans(const scalar p, const scalar Tt) const;

        //- Thermal conductivity, rotational mode [W/mK]
        inline scalar kappar(const scalar p, const scalar Tr) const;

        //- Thermal conductivity, trans-rotational mode [W/mK]
        inline scalar kappatr(const scalar p, const scalar Tt) const;

        //- Thermal conductivity, vibro-electronic mode [W/mK]
        inline scalar kappave
        (
            const scalar p,
            const scalar Tt,
            const scalar Tve
        ) const;

        //- Thermal diffusivity, trans-rotational mode [kg/ms]
        inline scalar alphatr(const scalar p, const scalar Tt) const;

        //- Thermal diffusivity, vib-electronic mode [kg/ms]
        inline scalar alphave
        (
            const scalar p,
            const scalar Tt,
            const scalar Tve
        ) const;

        //- Write to Ostream
        void write(Ostream& os) const;


    // Member operators

        inline constantTransport& operator=(const constantTransport&);

        inline void operator+=(const constantTransport&);

        inline void operator-=(const constantTransport&);

        inline void operator*=(const scalar);


    // Friend operators

        friend constantTransport operator+ <Thermo>
        (
            const constantTransport&,
            const constantTransport&
        );

        friend constantTransport operator- <Thermo>
        (
            const constantTransport&,
            const constantTransport&
        );

        friend constantTransport operator* <Thermo>
        (
            const scalar,
            const constantTransport&
        );

        friend constantTransport operator== <Thermo>
        (
            const constantTransport&,
            const constantTransport&
        );


    // Ostream Operator

        friend Ostream& operator<< <Thermo>
        (
            Ostream&,
            const constantTransport&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "constantTransportI.H"

#ifdef NoRepository
#   include "constantTransport.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
