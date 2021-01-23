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
    Foam::multiThermo

Description
    Basic thermodynamics type based on the use of fitting functions for
    cp, h, s obtained from the template argument type multiThermo.  All other
    properties are derived from these primitive functions.

SourceFiles
    multiThermoI.H
    multiThermo.C

\*---------------------------------------------------------------------------*/

#ifndef multiThermo_H
#define multiThermo_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace species
{

// Forward declaration of friend functions and operators

template<class MultiThermo, template<class> class Type> class multiThermo;

template<class MultiThermo, template<class> class Type>
inline multiThermo<MultiThermo, Type> operator+
(
    const multiThermo<MultiThermo, Type>&,
    const multiThermo<MultiThermo, Type>&
);

template<class MultiThermo, template<class> class Type>
inline multiThermo<MultiThermo, Type> operator-
(
    const multiThermo<MultiThermo, Type>&,
    const multiThermo<MultiThermo, Type>&
);

template<class MultiThermo, template<class> class Type>
inline multiThermo<MultiThermo, Type> operator*
(
    const scalar,
    const multiThermo<MultiThermo, Type>&
);

template<class MultiThermo, template<class> class Type>
inline multiThermo<MultiThermo, Type> operator==
(
    const multiThermo<MultiThermo, Type>&,
    const multiThermo<MultiThermo, Type>&
);

template<class MultiThermo, template<class> class Type>
Ostream& operator<<
(
    Ostream&,
    const multiThermo<MultiThermo, Type>&
);


/*---------------------------------------------------------------------------*\
                        Class multiThermo Declaration
\*---------------------------------------------------------------------------*/

template<class MultiThermo, template<class> class Type>
class multiThermo
:
    public MultiThermo,
    public Type<multiThermo<MultiThermo, Type> >
{
    // Private data

        //- Convergence tolerance of energy -> temperature inversion functions
        static const scalar tol_;

        //- Max number of iterations in energy->temperature inversion functions
        static const int maxIter_;


    // Private Member Functions

        //- Return the temperature corresponding to the value of the
        //  thermodynamic property f, given the function f = F(p, T)
        //  and dF(p, T)/dT
        //  Root-finding algorithm : Newton's method
        inline scalar T
        (
            scalar f,
            scalar p,
            scalar T0,
            scalar (multiThermo::*F)(const scalar, const scalar) const,
            scalar (multiThermo::*dFdT)(const scalar, const scalar) const,
            scalar (multiThermo::*limit)(const scalar) const
        ) const;

        // Function overload for mode vibMode
        inline scalar T
        (
            label vibMode,
            scalar f,
            scalar p,
            scalar T0,
            scalar (multiThermo::*F)
                (
                    const label,
                    const scalar,
                    const scalar
                ) const,
            scalar (multiThermo::*dFdT)
                (
                    const label,
                    const scalar,
                    const scalar
                ) const,
            scalar (multiThermo::*limit)(const scalar) const
        ) const;


public:

    // Constructors

        //- construct from components
        inline multiThermo(const MultiThermo& sp);

        //- Construct from Istream
        multiThermo(Istream&);

        //- Construct from dictionary
        multiThermo(const dictionary& dict);

        //- Construct as named copy
        inline multiThermo(const word& name, const multiThermo&);


    // Member Functions

        //- Return the instantiated type name
        static word typeName()
        {
            return
                  MultiThermo::typeName() + ','
                + Type<multiThermo<MultiThermo, Type> >::typeName();
        }

        // Fundamental properties
        // (These functions must be provided in derived types)

            // Heat capacity at constant volume [J/(kmol K)]
            // trans-rotational energy mode
            // scalar cv_t(const scalar, const scalar) const;

            // Heat capacity at constant volume [J/(kmol K)]
            // rotational energy mode
            // scalar cv_r(const scalar, const scalar) const;

            // Heat capacity at constant volume [J/(kmol K)]
            // vibrational energy mode, mode vibMode
            // scalar cv_v_mode(const label, const scalar, const scalar) const;

            // Heat capacity at constant volume [J/(kmol K)]
            // vibrational energy mode
            // scalar cv_v(const scalar, const scalar) const;

            // Heat capacity at constant volume [J/(kmol K)]
            // electronic energy mode
            // scalar cv_el(const scalar, const scalar) const;

            // Heat capacity at constant volume [J/(kmol K)]
            // electron energy mode
            // scalar cv_e(const scalar, const scalar) const;

            // Absolute enthalpy [J/kmol]
            // scalar ha(const scalar) const;

            // Sensible internal energy [J/kmol], trans-rotational energy mode
            // scalar ets(const scalar, const scalar) const;

            // Sensible internal energy [J/kmol], rotational energy mode
            // scalar ers(const scalar, const scalar) const;

            // Sensible internal energy [J/kmol]
            // for a given vibrational energy mode
            // scalar evs_mode(const label, const scalar, const scalar) const;

            // Sensible internal energy [J/kmol], vibrational energy mode
            // scalar evs(const scalar, const scalar) const;

            // Sensible internal energy [J/kmol], electronic energy mode
            // scalar eels(const scalar, const scalar) const;

            // Sensible internal energy [J/kmol], electron mode
            // scalar ees(const scalar, const scalar) const;

            //- Degrees of freedom associated with the rotational energy mode
            // inline scalar zetar
            // (
            //    const scalar,
            //    const scalar,
            //    const scalar
            // ) const;

            //- Degrees of freedom associated with the vibrational energy mode, mode vibMode
            // inline scalar zetav_mode(const label, const scalar, const scalar) const;

            //- Degrees of freedom associated with the vibrational energy mode
            // inline scalar zetav(const scalar, const scalar) const;

            //- Degrees of freedom associated with the electronic energy mode
            // inline scalar zetael(const scalar, const scalar) const;

            //- Degrees of freedom associated with the electron energy mode
            // inline scalar zetae(const scalar, const scalar) const;

            // Chemical enthalpy [J/kmol]
            // scalar hc() const;

            // Entropy [J/(kmol K)]
            // scalar s(const scalar, const scalar, const scalar) const;


        // Calculate and return derived properties
        // (These functions do not to be provided in derived types)

            // Mole specific properties

                //- Name of enthalpy/internal energy
                static inline word heName();

                //- Enthalpy/Internal energy [J/kmol] (backward compatibility)
                inline scalar he(const scalar p, const scalar T) const;
                
                //- Enthalpy/Internal energy [J/kmol]
                inline scalar he
                (
                    const scalar p,
                    const scalar T,
                    const scalar Tv
                ) const;

                //- Trans-rotational enthalpy/internal energy [J/kmol]
                inline scalar het(const scalar p, const scalar T) const;

                //- Vibrational enthalpy/internal energy [J/kmol]
                inline scalar hev(const scalar p, const scalar Tv) const;

                //- Vibrational enthalpy/internal energy [J/kmol], mode vibMode
                inline scalar hev_mode
                (
                    const label vibMode,
                    const scalar p,
                    const scalar Tv
                ) const;

                //- Electronic enthalpy/internal energy [J/kmol]
                inline scalar heel(const scalar p, const scalar Tv) const;

                //- Vibro-electronic enthalpy/internal energy [J/kmol]
                inline scalar hevel(const scalar p, const scalar Tv) const;

                //- Heat capacity at constant volume [J/(kmol K)]
                //  (backward compatibility)
                inline scalar cv(const scalar p, const scalar T) const;
                
                //- Heat capacity at constant volume [J/(kmol K)]
                inline scalar cv
                (
                    const scalar p,
                    const scalar T,
                    const scalar Tv
                ) const;

                //- Heat capacity at constant volume [J/(kmol K)]
                //  vibro-electronic energy mode
                inline scalar cv_vel(const scalar p, const scalar Tv) const;

                //- Heat capacity at constant pressure [J/(kmol K)]
                //  (backward compatibility)
                inline scalar cp(const scalar p, const scalar T) const;
                
                //- Heat capacity at constant pressure [J/(kmol K)]
                inline scalar cp
                (
                    const scalar p,
                    const scalar T,
                    const scalar Tv
                ) const;

                //- Heat capacity at constant pressure [J/(kmol K)]
                //  trans-rotational energy mode
                inline scalar cp_t(const scalar p, const scalar T) const;

                //- Heat capacity at constant pressure [J/(kmol K)]
                //  vibrational energy mode vibMode 
                inline scalar cp_v_mode
                (
                    const label vibMode,
                    const scalar p,
                    const scalar Tv
                ) const;

                //- Heat capacity at constant pressure [J/(kmol K)]
                //  vibrational energy mode 
                inline scalar cp_v(const scalar p, const scalar Tv) const;

                //- Heat capacity at constant pressure [J/(kmol K)]
                //  electronic energy mode
                inline scalar cp_el(const scalar p, const scalar Tv) const;

                //- gamma = cp/cv []
                inline scalar gamma
                (
                    const scalar p,
                    const scalar T,
                    const scalar Tv = 0
                ) const;

                //- Sensible internal energy [J/kmol] (backward compatibility)
                inline scalar es(const scalar p, const scalar T) const;
                
                //- Sensible internal energy [J/kmol]
                inline scalar es
                (
                    const scalar p,
                    const scalar T,
                    const scalar Tv
                ) const;

                //- Sensible enthalpy [J/kmol]
                //  vibro-electronic energy mode vibMode
                inline scalar evels_mode
                (
                    const label vibMode,
                    const scalar p,
                    const scalar Tv
                ) const;

                //- Sensible enthalpy [J/kmol], vibro-electronic energy mode
                inline scalar evels(const scalar p, const scalar Tv) const;

                //- Sensible enthalpy [J/kmol]
                inline scalar hs(const scalar p, const scalar T) const;
                inline scalar hs
                (
                    const scalar p,
                    const scalar T,
                    const scalar Tv
                ) const;

                //- Sensible enthalpy [J/kmol], trans-rotational energy mode
                inline scalar hts(const scalar p, const scalar T) const;

                //- Sensible enthalpy [J/kmol], vibrational energy mode vibMode
                inline scalar hvs_mode
                (
                    const label vibMode,
                    const scalar p,
                    const scalar Tv
                ) const;

                //- Sensible enthalpy [J/kmol], vibrational energy mode
                inline scalar hvs(const scalar p, const scalar Tv) const;

                //- Sensible enthalpy [J/kmol], electronic energy mode
                inline scalar hels(const scalar p, const scalar Tv) const;

                //- Absolute internal energy [J/kmol]
                inline scalar ea
                (
                    const scalar p,
                    const scalar T,
                    const scalar Tv = 0
                ) const;

                //- Gibbs free energy [J/kmol]
                inline scalar g
                (
                    const scalar p,
                    const scalar T,
                    const scalar Tv = 0
                ) const;

                //- Helmholtz free energy [J/kmol]
                inline scalar a
                (
                    const scalar p,
                    const scalar T,
                    const scalar Tv = 0
                ) const;


            // Mass specific properties

                //- Heat capacity at constant pressure [J/(kg K)]
                inline scalar Cp
                (
                    const scalar p,
                    const scalar T,
                    const scalar Tv = 0
                ) const;

                //- Heat capacity at constant volume [J/(kg K)]
                //  (backward compatibility)
                inline scalar Cv(const scalar p, const scalar T) const;
                
                //- Heat capacity at constant volume [J/(kg K)]
                inline scalar Cv
                (
                    const scalar p,
                    const scalar T,
                    const scalar Tv
                ) const;

                //- Heat capacity at constant pressure/volume [J/(kg K)]
                inline scalar Cpv
                (
                    const scalar p,
                    const scalar T,
                    const scalar Tv = 0
                ) const;

                //- Heat capacity at constant volume [J/(kg K)]
                //  translational energy mode
                inline scalar Cv_trans(const scalar p, const scalar T) const;

                //- Heat capacity at constant volume [J/(kg K)]
                //  rotational energy mode
                inline scalar Cv_r(const scalar p, const scalar T) const;

                //- Heat capacity at constant volume for the
                //  trans-rotational energy mode [J/(kg K)]
                inline scalar Cv_t(const scalar p, const scalar T) const;

                //- Heat capacity at constant volume [J/(kg K)]
                //  vibrational energy mode vibMode
                inline scalar Cv_v_mode
                (
                    const label vibMode,
                    const scalar p,
                    const scalar Tv
                ) const;

                //- Heat capacity at constant volume [J/(kg K)]
                //  vibrational energy mode
                inline scalar Cv_v(const scalar p, const scalar Tv) const;

                //- Heat capacity at constant volume [J/(kg K)]
                //  electronic energy mode
                inline scalar Cv_el(const scalar p, const scalar Tv) const;

                //- Heat capacity at constant volume [J/(kg K)]
                //  vibro-electronic energy mode
                inline scalar Cv_vel(const scalar p, const scalar Tv) const;

                //- Heat capacity at constant pressure [J/(kg K)]
                //  trans-rotational energy mode
                inline scalar Cp_t(const scalar p, const scalar T) const;

                //- Heat capacity at constant pressure [J/(kg K)]
                //  vibrational energy mode
                inline scalar Cp_v_mode
                (
                    const label vibMode,
                    const scalar p,
                    const scalar Tv
                ) const;

                //- Heat capacity at constant pressure [J/(kg K)]
                //  vibrational energy mode
                inline scalar Cp_v(const scalar p, const scalar Tv) const;

                //- Heat capacity at constant pressure [J/(kg K)]
                //  electronic energy mode
                inline scalar Cp_el(const scalar p, const scalar Tv) const;

                //- Enthalpy/Internal energy [J/kg] (backward compatibility)
                inline scalar HE(const scalar p, const scalar T) const;
                
                //- Enthalpy/Internal energy [J/kg]
                inline scalar HE
                (
                    const scalar p,
                    const scalar T,
                    const scalar Tv
                ) const;

                //- Trans-rotational enthalpy/internal energy [J/kg]
                inline scalar HEt(const scalar p, const scalar T) const;

                //- Vibrational enthalpy/internal energy [J/kg]
                inline scalar HEv(const scalar p, const scalar Tv) const;

                //- Vibrational enthalpy/internal energy [J/kg], mode vibMode
                inline scalar HEv_mode
                (
                    const label vibMode,
                    const scalar p,
                    const scalar Tv
                ) const;

                //- Electronic enthalpy/internal energy [J/kg]
                inline scalar HEel(const scalar p, const scalar Tv) const;

                //- Vibro-electronic enthalpy/internal energy [J/kg]
                inline scalar HEvel(const scalar p, const scalar Tv) const;

                //- Sensible enthalpy [J/kg]
                inline scalar Hs
                (
                    const scalar p,
                    const scalar T,
                    const scalar Tv = 0
                ) const;

                //- Sensible enthalpy [J/kg], trans-rotational energy mode
                inline scalar Hts(const scalar p, const scalar T) const;

                //- Sensible enthalpy [J/kg], vibrational energy mode vibMode
                inline scalar Hvs_mode
                (
                    const label vibMode,
                    const scalar p,
                    const scalar Tv
                ) const;

                //- Sensible enthalpy [J/kg], vibrational energy mode
                inline scalar Hvs(const scalar p, const scalar Tv) const;

                //- Sensible enthalpy [J/kg], electronic energy mode
                inline scalar Hels(const scalar p, const scalar Tv) const;

                //- Chemical enthalpy [J/kg]
                inline scalar Hc() const;

                //- Absolute Enthalpy [J/kg]
                inline scalar Ha
                (
                    const scalar p,
                    const scalar T,
                    const scalar Tv = 0
                ) const;

                //- Entropy [J/(kg K)]
                inline scalar S
                (
                    const scalar p,
                    const scalar T,
                    const scalar Tv = 0
                ) const;

                //- Internal energy [J/kg]
                inline scalar E
                (
                    const scalar p,
                    const scalar T,
                    const scalar Tv = 0
                ) const;

                //- Sensible internal energy [J/kg] (backward compatibility)
                inline scalar Es(const scalar p, const scalar T) const;
                
                //- Sensible internal energy [J/kg]
                inline scalar Es
                (
                    const scalar p,
                    const scalar T,
                    const scalar Tv
                ) const;

                //- Sensible internal energy [J/kg]
                //  trans-rotational energy mode
                inline scalar Ets(const scalar p, const scalar T) const;

                //- Sensible internal energy [J/kg]
                //  vibrational energy mode vibMode
                inline scalar Evs_mode
                (
                    const label vibMode,
                    const scalar p,
                    const scalar Tv
                ) const;

                //- Sensible internal energy [J/kg], vibrational energy mode
                inline scalar Evs(const scalar p, const scalar Tv) const;

                //- Sensible internal energy [J/kg], electronic energy mode
                inline scalar Eels(const scalar p, const scalar Tv) const;

                //- Sensible internal energy [J/kg]
                //  vibro-electronic energy mode
                inline scalar Evels(const scalar p, const scalar Tv) const;

                //- Absolute internal energy [J/kg]
                inline scalar Ea
                (
                    const scalar p,
                    const scalar T,
                    const scalar Tv = 0
                ) const;

                //- Gibbs free energy [J/kg]
                inline scalar G
                (
                    const scalar p,
                    const scalar T,
                    const scalar Tv = 0
                ) const;

                //- Helmholtz free energy [J/kg]
                inline scalar A
                (
                    const scalar p,
                    const scalar T,
                    const scalar Tv = 0
                ) const;


        // Equilibrium reaction thermodynamics

            //- Equilibrium constant [] i.t.o fugacities
            //  = PIi(fi/Pstd)^nui
            inline scalar K(const scalar p, const scalar T) const;

            //- Equilibrium constant [] i.t.o. partial pressures
            //  = PIi(pi/Pstd)^nui
            //  For low pressures (where the gas mixture is near perfect) Kp = K
            inline scalar Kp(const scalar p, const scalar T) const;

            //- Equilibrium constant i.t.o. molar concentration
            //  = PIi(ci/cstd)^nui
            //  For low pressures (where the gas mixture is near perfect)
            //  Kc = Kp(pstd/(RR*T))^nu
            inline scalar Kc(const scalar p, const scalar T) const;

            //- Equilibrium constant [] i.t.o. mole-fractions
            //  For low pressures (where the gas mixture is near perfect)
            //  Kx = Kp(pstd/p)^nui
            inline scalar Kx
            (
                const scalar p,
                const scalar T
            ) const;

            //- Equilibrium constant [] i.t.o. number of moles
            //  For low pressures (where the gas mixture is near perfect)
            //  Kn = Kp(n*pstd/p)^nui where n = number of moles in mixture
            inline scalar Kn
            (
                const scalar p,
                const scalar T,
                const scalar n
            ) const;


        // Energy->temperature  inversion functions

            //- Trans-rotational temperature from enthalpy or internal
            //  energy given an initial temperature T0
            inline scalar TtHEt
            (
                const scalar H,
                const scalar p,
                const scalar T0
            ) const;

            //- Vibrational temperature from enthalpy or internal energy
            //  given an initial temperature T0
            inline scalar TvHEv
            (
                const scalar H,
                const scalar p,
                const scalar T0
            ) const;

            //- Vibrational temperature from enthalpy or internal energy
            //  given an initial temperature T0, mode vibMode
            inline scalar TvHEv_mode
            (
                const label vibMode,
                const scalar H,
                const scalar p,
                const scalar T0
            ) const;

            //- Vibro-electronic temperature from enthalpy or internal energy
            //  given an initial temperature T0
            inline scalar TvelHEvel
            (
                const scalar H,
                const scalar p,
                const scalar T0
            ) const;

            //- Trans-rotational temperature from sensible enthalpy
            //  given an initial T0
            inline scalar TtHts
            (
                const scalar Hs,
                const scalar p,
                const scalar T0
            ) const;

            //- Vibrational temperature from sensible enthalpy given an
            //  initial T0
            inline scalar TvHvs
            (
                const scalar Hs,
                const scalar p,
                const scalar T0
            ) const;

            //- Temperature from sensible internal energy
            //  given an initial temperature T0
            inline scalar TEs
            (
                const scalar E,
                const scalar p,
                const scalar T0
            ) const;

            //- Trans-rotational temperature from sensible internal energy
            //  given an initial T0
            inline scalar TtEts
            (
                const scalar Es,
                const scalar p,
                const scalar T0
            ) const;

            //- Vibrational temperature from sensible internal energy given an
            //  initial T0, mode vibMode
            inline scalar TvEvs_mode
            (
                const label vibMode,
                const scalar Es_mode,
                const scalar p,
                const scalar T0
            ) const;

            //- Vibrational temperature from sensible internal energy given an
            //  initial T0
            inline scalar TvEvs
            (
                const scalar Es,
                const scalar p,
                const scalar T0
            ) const;

            //- Vibro-electronic temperature from sensible internal energy
            //  given an initial T0
            inline scalar TvelEvels
            (
                const scalar Es,
                const scalar p,
                const scalar T0
            ) const;


        // I-O

            //- Write to Ostream
            void write(Ostream& os) const;


    // Member operators

        inline void operator+=(const multiThermo&);
        inline void operator-=(const multiThermo&);

        inline void operator*=(const scalar);


    // Friend operators

        friend multiThermo operator+ <MultiThermo, Type>
        (
            const multiThermo&,
            const multiThermo&
        );

        friend multiThermo operator- <MultiThermo, Type>
        (
            const multiThermo&,
            const multiThermo&
        );

        friend multiThermo operator* <MultiThermo, Type>
        (
            const scalar s,
            const multiThermo&
        );

        friend multiThermo operator== <MultiThermo, Type>
        (
            const multiThermo&,
            const multiThermo&
        );


    // Ostream Operator

        friend Ostream& operator<< <MultiThermo, Type>
        (
            Ostream&,
            const multiThermo&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace species
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "multiThermoI.H"

#ifdef NoRepository
#   include "multiThermo.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
