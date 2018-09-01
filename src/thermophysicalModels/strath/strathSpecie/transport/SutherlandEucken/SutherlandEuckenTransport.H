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
    Foam::SutherlandEuckenTransport

Description
    Transport package using Sutherland's formula.

    Templated into a given thermodynamics package (needed for thermal
    conductivity).

    Dynamic viscosity [kg/m.s]
    \f[
        \mu = A_s \frac{\sqrt{T}}{1 + T_s / T}
    \f]

SourceFiles
    SutherlandEuckenTransportI.H
    SutherlandEuckenTransport.C

\*---------------------------------------------------------------------------*/

#ifndef SutherlandEuckenTransport_H
#define SutherlandEuckenTransport_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

template<class Thermo> class SutherlandEuckenTransport;

template<class Thermo>
inline SutherlandEuckenTransport<Thermo> operator+
(
    const SutherlandEuckenTransport<Thermo>&,
    const SutherlandEuckenTransport<Thermo>&
);

template<class Thermo>
inline SutherlandEuckenTransport<Thermo> operator-
(
    const SutherlandEuckenTransport<Thermo>&,
    const SutherlandEuckenTransport<Thermo>&
);

template<class Thermo>
inline SutherlandEuckenTransport<Thermo> operator*
(
    const scalar,
    const SutherlandEuckenTransport<Thermo>&
);

template<class Thermo>
inline SutherlandEuckenTransport<Thermo> operator==
(
    const SutherlandEuckenTransport<Thermo>&,
    const SutherlandEuckenTransport<Thermo>&
);

template<class Thermo>
Ostream& operator<<
(
    Ostream&,
    const SutherlandEuckenTransport<Thermo>&
);


/*---------------------------------------------------------------------------*\
                    Class SutherlandEuckenTransport Declaration
\*---------------------------------------------------------------------------*/

template<class Thermo>
class SutherlandEuckenTransport
:
    public Thermo
{
    // Private data

        //- Sutherland's coefficients
        scalar As_, Ts_;
        
        //- Factor in vibrational thermal conductivity kappave
        // Hard-sphere: (7-2*omega_m)/5 with omega_m = 0.5
        // Default value: 1.2
        scalar eta_s_;


    // Private Member Functions

        //- Calculate the Sutherland coefficients
        //  given two viscosities and temperatures
        inline void calcCoeffs
        (
            const scalar mu1, const scalar T1,
            const scalar mu2, const scalar T2
            
        );


public:

    // Constructors

        //- Construct from components
        inline SutherlandEuckenTransport
        (
            const Thermo& t,
            const scalar As,
            const scalar Ts,
            const scalar eta_s
        );

        //- Construct from two viscosities
        inline SutherlandEuckenTransport
        (
            const Thermo& t,
            const scalar mu1, const scalar T1,
            const scalar mu2, const scalar T2
        );

        //- Construct as named copy
        inline SutherlandEuckenTransport(const word&, const SutherlandEuckenTransport&);

        //- Construct from Istream
        SutherlandEuckenTransport(Istream&);

        //- Construct from dictionary
        SutherlandEuckenTransport(const dictionary& dict);

        //- Construct and return a clone
        inline autoPtr<SutherlandEuckenTransport> clone() const;

        // Selector from Istream
        inline static autoPtr<SutherlandEuckenTransport> New(Istream& is);

        // Selector from dictionary
        inline static autoPtr<SutherlandEuckenTransport> New(const dictionary& dict);


    // Member functions

        //- Return the instantiated type name
        static word typeName()
        {
            return "SutherlandEucken<" + Thermo::typeName() + '>';
        }

        //- Dynamic viscosity [kg/ms]
        inline scalar mu(const scalar p, const scalar T) const;

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

        inline SutherlandEuckenTransport& operator=(const SutherlandEuckenTransport&);

        inline void operator+=(const SutherlandEuckenTransport&);

        inline void operator-=(const SutherlandEuckenTransport&);

        inline void operator*=(const scalar);


    // Friend operators

        friend SutherlandEuckenTransport operator+ <Thermo>
        (
            const SutherlandEuckenTransport&,
            const SutherlandEuckenTransport&
        );

        friend SutherlandEuckenTransport operator- <Thermo>
        (
            const SutherlandEuckenTransport&,
            const SutherlandEuckenTransport&
        );

        friend SutherlandEuckenTransport operator* <Thermo>
        (
            const scalar,
            const SutherlandEuckenTransport&
        );

        friend SutherlandEuckenTransport operator== <Thermo>
        (
            const SutherlandEuckenTransport&,
            const SutherlandEuckenTransport&
        );


    // Ostream Operator

        friend Ostream& operator<< <Thermo>
        (
            Ostream&,
            const SutherlandEuckenTransport&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "SutherlandEuckenTransportI.H"

#ifdef NoRepository
#   include "SutherlandEuckenTransport.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
