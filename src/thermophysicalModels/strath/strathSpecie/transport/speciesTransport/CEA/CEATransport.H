/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016-2021 hyStrath
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
    Foam::CEATransport

Description
    Transport package using CEA formulas.

    Templated into a given thermodynamics package (needed for thermal
    conductivity).

    Dynamic viscosity [kg/m.s]
    \f[
        \mu_k = exp(Amu_k*ln(T) + Bmu_k/T + Cmu_k/T^2 + Dmu_k)
    \f]

    Thermal conductivity [W/mK]
    \f[
        \kappa_tr_k = exp(Akp_k*ln(T) + Bkp_k/T + Ckp_k/T^2 + Dkp_k)
    \f]

    \kappa_ve_k amended to be eta_s_k*mu_k*Cv_ve_k

SourceFiles
    CEATransportI.H
    CEATransport.C

\*---------------------------------------------------------------------------*/

#ifndef CEATransport_H
#define CEATransport_H

#include "DynamicList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

template<class Thermo> class CEATransport;

template<class Thermo>
inline CEATransport<Thermo> operator+
(
    const CEATransport<Thermo>&,
    const CEATransport<Thermo>&
);

template<class Thermo>
inline CEATransport<Thermo> operator-
(
    const CEATransport<Thermo>&,
    const CEATransport<Thermo>&
);

template<class Thermo>
inline CEATransport<Thermo> operator*
(
    const scalar,
    const CEATransport<Thermo>&
);

template<class Thermo>
inline CEATransport<Thermo> operator==
(
    const CEATransport<Thermo>&,
    const CEATransport<Thermo>&
);

template<class Thermo>
Ostream& operator<<
(
    Ostream&,
    const CEATransport<Thermo>&
);


/*---------------------------------------------------------------------------*\
                    Class CEATransport Declaration
\*---------------------------------------------------------------------------*/

template<class Thermo>
class CEATransport
:
    public Thermo
{

public:

    static const int nCoeffs_ = 4;
    typedef FixedList<scalar, nCoeffs_> CEATransportArray;

private:

    // Private data

        //- Temperature bounds for the viscosity and thermal conductivity lists
        DynamicList<scalar> temp_;

        //- CEA coefficients for viscosity and thermal conductivity
        DynamicList<CEATransportArray> mu_, kappa_;

        //- Factor in vibrational thermal conductivity kappave
        // Hard-sphere: (7-2*omega_m)/5 with omega_m = 0.5
        // Default value: 1.2
        scalar eta_s_;


public:

    // Constructors

        //- Construct from components
        inline CEATransport
        (
            const Thermo& t,
            const DynamicList<scalar>& temp,
            const DynamicList<CEATransportArray>& mu,
            const DynamicList<CEATransportArray>& kappa,
            const scalar eta_s
        );

        //- Construct as named copy
        inline CEATransport(const word&, const CEATransport&);

        //- Construct from Istream
        CEATransport(Istream&);

        //- Construct from dictionary
        CEATransport(const dictionary& dict);

        //- Construct and return a clone
        inline autoPtr<CEATransport> clone() const;

        // Selector from Istream
        inline static autoPtr<CEATransport> New(Istream& is);

        // Selector from dictionary
        inline static autoPtr<CEATransport> New(const dictionary& dict);


    // Member functions

        //- Return the instantiated type name
        static word typeName()
        {
            return "CEA<" + Thermo::typeName() + '>';
        }

        //- Return the viscosity coefficients for a given temperature
        inline const CEATransportArray& viscosityCoeffs(const scalar T) const;

        //- Return the thermal conductivity coefficients for a given temperature
        inline const CEATransportArray& thermalConductivityCoeffs
        (
            const scalar T
        ) const;

        //- Dynamic viscosity [kg/ms]
        inline scalar mu(const scalar p, const scalar Tt) const;

        //- Thermal conductivity, trans-rotational mode [W/mK]
        inline scalar kappatr(const scalar p, const scalar Tt) const;

        //- Thermal conductivity, vibro-electronic energy mode [W/mK]
        inline scalar kappave
        (
            const scalar p,
            const scalar Tt,
            const scalar Tve
        ) const;

        //- Write to Ostream
        void write(Ostream& os) const;


    // Member operators

        inline CEATransport& operator=(const CEATransport&);

        inline void operator+=(const CEATransport&);

        inline void operator-=(const CEATransport&);

        inline void operator*=(const scalar);


    // Friend operators

        friend CEATransport operator+ <Thermo>
        (
            const CEATransport&,
            const CEATransport&
        );

        friend CEATransport operator- <Thermo>
        (
            const CEATransport&,
            const CEATransport&
        );

        friend CEATransport operator* <Thermo>
        (
            const scalar,
            const CEATransport&
        );

        friend CEATransport operator== <Thermo>
        (
            const CEATransport&,
            const CEATransport&
        );


    // Ostream Operator

        friend Ostream& operator<< <Thermo>
        (
            Ostream&,
            const CEATransport&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "CEATransportI.H"

#ifdef NoRepository
#   include "CEATransport.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
