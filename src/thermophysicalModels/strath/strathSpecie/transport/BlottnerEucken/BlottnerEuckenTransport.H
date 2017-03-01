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
    Foam::BlottnerEuckenTransport

Description
    Transport package using Blottner and Eucken's formulas.

    Templated into a given thermodynamics package (needed for thermal
    conductivity).

    Dynamic viscosity [kg/m.s] - Blottner
    \f[
        \mu_k = 0.1*exp((A_k*ln(Tt)+B_k)*ln(Tt)+C_k)
    \f]
    
    Thermal conductivity [W/mK] - Eucken
    \f[
        \kappa_t_k = 3/2*mu_k*Cv_t_k
        \kappa_r_k = mu_k*Cv_r_k
        \kappa_ve_k = mu_k*Cv_ve_k
    \f]

SourceFiles
    BlottnerEuckenTransportI.H
    BlottnerEuckenTransport.C

\*---------------------------------------------------------------------------*/

#ifndef BlottnerEuckenTransport_H
#define BlottnerEuckenTransport_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

template<class Thermo> class BlottnerEuckenTransport;

template<class Thermo>
inline BlottnerEuckenTransport<Thermo> operator+
(
    const BlottnerEuckenTransport<Thermo>&,
    const BlottnerEuckenTransport<Thermo>&
);

template<class Thermo>
inline BlottnerEuckenTransport<Thermo> operator-
(
    const BlottnerEuckenTransport<Thermo>&,
    const BlottnerEuckenTransport<Thermo>&
);

template<class Thermo>
inline BlottnerEuckenTransport<Thermo> operator*
(
    const scalar,
    const BlottnerEuckenTransport<Thermo>&
);

template<class Thermo>
inline BlottnerEuckenTransport<Thermo> operator==
(
    const BlottnerEuckenTransport<Thermo>&,
    const BlottnerEuckenTransport<Thermo>&
);

template<class Thermo>
Ostream& operator<<
(
    Ostream&,
    const BlottnerEuckenTransport<Thermo>&
);


/*---------------------------------------------------------------------------*\
                    Class BlottnerEuckenTransport Declaration
\*---------------------------------------------------------------------------*/

template<class Thermo>
class BlottnerEuckenTransport
:
    public Thermo
{
    // Private data

        // Blottner's coefficients
        scalar Ak_, Bk_, Ck_;


public:

    // Constructors

        //- Construct from components
        inline BlottnerEuckenTransport
        (
            const Thermo& t,
            const scalar Ak,
            const scalar Bk,
            const scalar Ck
        );

        //- Construct as named copy
        inline BlottnerEuckenTransport(const word&, const BlottnerEuckenTransport&);

        //- Construct from Istream
        BlottnerEuckenTransport(Istream&);

        //- Construct from dictionary
        BlottnerEuckenTransport(const dictionary& dict);

        //- Construct and return a clone
        inline autoPtr<BlottnerEuckenTransport> clone() const;

        // Selector from Istream
        inline static autoPtr<BlottnerEuckenTransport> New(Istream& is);

        // Selector from dictionary
        inline static autoPtr<BlottnerEuckenTransport> New(const dictionary& dict);


    // Member functions

        //- Return the instantiated type name
        static word typeName()
        {
            return "BlottnerEucken<" + Thermo::typeName() + '>';
        }

        //- Dynamic viscosity [kg/ms]
        inline scalar mu(const scalar p, const scalar Tt) const;

        //- Thermal conductivity, translational mode [W/mK]
        inline scalar kappatrans(const scalar p, const scalar Tt) const;
        
        //- Thermal conductivity, rotational mode [W/mK]
        inline scalar kappar(const scalar p, const scalar Tr) const;
        
        //- Thermal conductivity, trans-rotational mode [W/mK]
        inline scalar kappatr(const scalar p, const scalar Tt) const;
        
        //- Thermal conductivity, vib-electronic mode [W/mK]
        inline scalar kappave(const scalar p, const scalar Tt, const scalar Tve) const;
        
        //- Thermal diffusivity, trans-rotational mode [kg/ms]
        inline scalar alphatr(const scalar p, const scalar Tt) const;
        
        //- Thermal diffusivity, vib-electronic mode [kg/ms]
        inline scalar alphave(const scalar p, const scalar Tt, const scalar Tve) const;

        //- Write to Ostream
        void write(Ostream& os) const;


    // Member operators

        inline BlottnerEuckenTransport& operator=(const BlottnerEuckenTransport&);

        inline void operator+=(const BlottnerEuckenTransport&);

        inline void operator-=(const BlottnerEuckenTransport&);

        inline void operator*=(const scalar);


    // Friend operators

        friend BlottnerEuckenTransport operator+ <Thermo>
        (
            const BlottnerEuckenTransport&,
            const BlottnerEuckenTransport&
        );

        friend BlottnerEuckenTransport operator- <Thermo>
        (
            const BlottnerEuckenTransport&,
            const BlottnerEuckenTransport&
        );

        friend BlottnerEuckenTransport operator* <Thermo>
        (
            const scalar,
            const BlottnerEuckenTransport&
        );

        friend BlottnerEuckenTransport operator== <Thermo>
        (
            const BlottnerEuckenTransport&,
            const BlottnerEuckenTransport&
        );


    // Ostream Operator

        friend Ostream& operator<< <Thermo>
        (
            Ostream&,
            const BlottnerEuckenTransport&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "BlottnerEuckenTransportI.H"

#ifdef NoRepository
#   include "BlottnerEuckenTransport.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
