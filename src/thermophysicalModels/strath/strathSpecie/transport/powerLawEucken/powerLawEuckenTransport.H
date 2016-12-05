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
    Foam::powerLawEuckenTransport

Description
    Transport package using Blottner and Eucken's formulas.

    Templated into a given thermodynamics package (needed for thermal
    conductivity).

    Dynamic viscosity [kg/m.s] - powerLaw
    \f[
        \mu_k = mu_ref_k*(Tt/T_ref)^omega_k
    \f]
    
    Thermal conductivity [W/m.K] - Eucken
    \f[
        \kappa_tr_k = 5/2*mu_k*Cv_tr_k
        \kappa_ve_k = mu_k*Cv_ve_k
    \f]

SourceFiles
    powerLawEuckenTransportI.H
    powerLawEuckenTransport.C

\*---------------------------------------------------------------------------*/

#ifndef powerLawEuckenTransport_H
#define powerLawEuckenTransport_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

template<class Thermo> class powerLawEuckenTransport;

template<class Thermo>
inline powerLawEuckenTransport<Thermo> operator+
(
    const powerLawEuckenTransport<Thermo>&,
    const powerLawEuckenTransport<Thermo>&
);

template<class Thermo>
inline powerLawEuckenTransport<Thermo> operator-
(
    const powerLawEuckenTransport<Thermo>&,
    const powerLawEuckenTransport<Thermo>&
);

template<class Thermo>
inline powerLawEuckenTransport<Thermo> operator*
(
    const scalar,
    const powerLawEuckenTransport<Thermo>&
);

template<class Thermo>
inline powerLawEuckenTransport<Thermo> operator==
(
    const powerLawEuckenTransport<Thermo>&,
    const powerLawEuckenTransport<Thermo>&
);

template<class Thermo>
Ostream& operator<<
(
    Ostream&,
    const powerLawEuckenTransport<Thermo>&
);


/*---------------------------------------------------------------------------*\
                    Class powerLawEuckenTransport Declaration
\*---------------------------------------------------------------------------*/

template<class Thermo>
class powerLawEuckenTransport
:
    public Thermo
{
    // Private data

        // Blottner's coefficients
        scalar dref_, omega_;


public:

    // Constructors

        //- Construct from components
        inline powerLawEuckenTransport
        (
            const Thermo& t,
            const scalar dref,
            const scalar omega
        );

        //- Construct as named copy
        inline powerLawEuckenTransport(const word&, const powerLawEuckenTransport&);

        //- Construct from Istream
        powerLawEuckenTransport(Istream&);

        //- Construct from dictionary
        powerLawEuckenTransport(const dictionary& dict);

        //- Construct and return a clone
        inline autoPtr<powerLawEuckenTransport> clone() const;

        // Selector from Istream
        inline static autoPtr<powerLawEuckenTransport> New(Istream& is);

        // Selector from dictionary
        inline static autoPtr<powerLawEuckenTransport> New(const dictionary& dict);


    // Member functions

        //- Return the instantiated type name
        static word typeName()
        {
            return "powerLawEucken<" + Thermo::typeName() + '>';
        }

        //- Dynamic viscosity [kg/ms]
        inline scalar mu(const scalar p, const scalar Tt) const;

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

        inline powerLawEuckenTransport& operator=(const powerLawEuckenTransport&);

        inline void operator+=(const powerLawEuckenTransport&);

        inline void operator-=(const powerLawEuckenTransport&);

        inline void operator*=(const scalar);


    // Friend operators

        friend powerLawEuckenTransport operator+ <Thermo>
        (
            const powerLawEuckenTransport&,
            const powerLawEuckenTransport&
        );

        friend powerLawEuckenTransport operator- <Thermo>
        (
            const powerLawEuckenTransport&,
            const powerLawEuckenTransport&
        );

        friend powerLawEuckenTransport operator* <Thermo>
        (
            const scalar,
            const powerLawEuckenTransport&
        );

        friend powerLawEuckenTransport operator== <Thermo>
        (
            const powerLawEuckenTransport&,
            const powerLawEuckenTransport&
        );


    // Ostream Operator

        friend Ostream& operator<< <Thermo>
        (
            Ostream&,
            const powerLawEuckenTransport&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "powerLawEuckenTransportI.H"

#ifdef NoRepository
#   include "powerLawEuckenTransport.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
