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
    Foam::decoupledEnergyModesThermo

Description
    DEM tables based thermodynamics package templated
    into the equation of state.

SourceFiles
    decoupledEnergyModesThermoI.H
    decoupledEnergyModesThermo.C

\*---------------------------------------------------------------------------*/

#ifndef decoupledEnergyModesThermo_H
#define decoupledEnergyModesThermo_H

#include "scalar.H"
#include "FixedList.H"
#include "DynamicList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

template<class EquationOfState> class decoupledEnergyModesThermo;

template<class EquationOfState>
inline decoupledEnergyModesThermo<EquationOfState> operator+
(
    const decoupledEnergyModesThermo<EquationOfState>&,
    const decoupledEnergyModesThermo<EquationOfState>&
);

template<class EquationOfState>
inline decoupledEnergyModesThermo<EquationOfState> operator-
(
    const decoupledEnergyModesThermo<EquationOfState>&,
    const decoupledEnergyModesThermo<EquationOfState>&
);

template<class EquationOfState>
inline decoupledEnergyModesThermo<EquationOfState> operator*
(
    const scalar,
    const decoupledEnergyModesThermo<EquationOfState>&
);

template<class EquationOfState>
inline decoupledEnergyModesThermo<EquationOfState> operator==
(
    const decoupledEnergyModesThermo<EquationOfState>&,
    const decoupledEnergyModesThermo<EquationOfState>&
);

template<class EquationOfState>
Ostream& operator<<
(
    Ostream&,
    const decoupledEnergyModesThermo<EquationOfState>&
);


/*---------------------------------------------------------------------------*\
                 Class decoupledEnergyModesThermo Declaration
\*---------------------------------------------------------------------------*/

template<class EquationOfState>
class decoupledEnergyModesThermo
:
    public EquationOfState
{

public:

    // Public data

        //- Number of coefficients in the input table
        static const int nCoeffs_ = 7;
        
        //- 
        // 0: cv translational
        // 1: cv rotational
        // 2: cv vibrational
        // 3: cv electronic
        // 4: cv electron
        // 5: chemical enthalpy
        // 6: chemical entropy
        typedef FixedList<scalar, nCoeffs_> cvCoeffsArray;
        
        //-
        typedef DynamicList<scalar> dynamicArray;


private:

    // Private data

        //- Temperature limits of applicability of functions
        scalar Tlow_, Thigh_;

        //-
        cvCoeffsArray decoupledCvCoeffs_;
        
        //-
        dynamicArray vibrationalList_, electronicList_;


    // Private Member Functions

        //- Check that input data is valid
        void checkInputData() const;

        //- Return the cv coefficients
        inline const cvCoeffsArray& cvCoeffs(const scalar T) const;

        //- Return the vibrational coefficients in the format
        //  (degeneracy - temperature)
        inline const dynamicArray& vibrationalCoeffs(const scalar T) const;

        //- Return the electronic coefficients in the format
        //  (degeneracy - temperature)
        inline const dynamicArray& electronicCoeffs(const scalar T) const;


public:

    // Constructors

        //- Construct from components
        inline decoupledEnergyModesThermo
        (
            const EquationOfState& st,
            const scalar Tlow,
            const scalar Thigh,
            const cvCoeffsArray& decoupledCvCoeffs,
            const dynamicArray& vibrationalList,
            const dynamicArray& electronicList
        );

        //- Construct from Istream
        decoupledEnergyModesThermo(Istream&);

        //- Construct from dictionary
        decoupledEnergyModesThermo(const dictionary& dict);

        //- Construct as a named copy
        inline decoupledEnergyModesThermo
        (
            const word&,
            const decoupledEnergyModesThermo&
        );


    // Member Functions

        //- Return the instantiated type name
        static word typeName()
        {
            return "decoupledEnergyModes<" + EquationOfState::typeName() + '>';
        }

        //- Limit the temperature to be in the range Tlow_ to Thigh_
        inline scalar limit(const scalar T) const;


        // Access

            //- Return const access to the low temperature limit
            inline scalar Tlow() const;

            //- Return const access to the high temperature limit
            inline scalar Thigh() const;

            //- Return const access to the list of cv coeffs
            inline const cvCoeffsArray& decoupledCvCoeffs() const;

            //- Return const access to the list of vibrational coeffs
            inline const dynamicArray& vibrationalList() const;

            //- Return const access to the list of electronic coeffs
            inline const dynamicArray& electronicList() const;



        // Fundamental properties

            //- Heat capacity at constant volume for the translational
            //  energy mode [J/(kmol K)]
            inline scalar cv_trans(const scalar p, const scalar T) const;

            //- Heat capacity at constant volume for the rotational
            //  energy mode [J/(kmol K)]
            inline scalar cv_r(const scalar p, const scalar T) const;

            //- Heat capacity at constant volume for the trans-rotational
            //  energy mode [J/(kmol K)]
            inline scalar cv_t(const scalar p, const scalar T) const;

            //- Heat capacity at constant volume for the vibrational
            //  energy mode, for the mode vibMode [J/(kmol K)]
            inline scalar cv_v_mode
            (
                const label vibMode,
                const scalar p,
                const scalar Tv
            ) const;

            //- Heat capacity at constant volume for the vibrational
            //  energy mode [J/(kmol K)]
            inline scalar cv_v(const scalar p, const scalar Tv) const;

            //- Heat capacity at constant volume for the electronic
            //  energy mode [J/(kmol K)]
            inline scalar cv_el(const scalar p, const scalar Tv) const;

            //- Heat capacity at constant volume for the electron energy
            //  mode [J/(kmol K)]
            inline scalar cv_e(const scalar p, const scalar Tv) const;

            //- Absolute enthalpy [J/kmol]
            inline scalar ha
            (
                const scalar p,
                const scalar T,
                const scalar Tv = 0
            ) const;

            //- Sensible trans-rotational energy [J/kmol]
            inline scalar ets(const scalar p, const scalar T) const;

            //- Sensible rotational energy [J/kmol]
            inline scalar ers(const scalar p, const scalar T) const;

            //- Sensible vibrational energy, for the vibrational mode
            //  vibMode [J/kmol]
            inline scalar evs_mode
            (
                const label vibMode,
                const scalar p,
                const scalar Tv
            ) const;

            //- Sensible vibrational energy [J/kmol]
            inline scalar evs(const scalar p, const scalar Tv) const;

            //- Sensible electronic energy [J/kmol]
            inline scalar eels(const scalar p, const scalar Tv) const;

             //- Sensible electron energy [J/kmol]
            inline scalar ees(const scalar p, const scalar Tv) const;

            //- Degrees of freedom associated with the rotational energy mode
            inline scalar zetar
            (
                const scalar p,
                const scalar T,
                const scalar Tv
            ) const;

            //- Degrees of freedom associated with the vibrational energy mode,
            //  for the vibrational mode vibMode
            inline scalar zetav_mode
            (
                const label vibMode,
                const scalar p,
                const scalar Tv
            ) const;

            //- Degrees of freedom associated with the vibrational energy mode
            inline scalar zetav(const scalar p, const scalar Tv) const;

            //- Degrees of freedom associated with the electronic energy mode
            inline scalar zetael(const scalar p, const scalar Tv) const;

            //- Degrees of freedom associated with the electron energy mode
            inline scalar zetae(const scalar p, const scalar Tv) const;

            //- Chemical enthalpy [J/kmol]
            inline scalar hc() const;

            //- Entropy [J/(kmol K)]
            inline scalar s
            (
                const scalar p,
                const scalar T,
                const scalar Tv = 0
            ) const;


        // I-O

            //- Write to Ostream
            void write(Ostream& os) const;


    // Member operators

        inline void operator+=(const decoupledEnergyModesThermo&);
        inline void operator-=(const decoupledEnergyModesThermo&);


    // Friend operators

        friend decoupledEnergyModesThermo operator+ <EquationOfState>
        (
            const decoupledEnergyModesThermo&,
            const decoupledEnergyModesThermo&
        );

        friend decoupledEnergyModesThermo operator- <EquationOfState>
        (
            const decoupledEnergyModesThermo&,
            const decoupledEnergyModesThermo&
        );

        friend decoupledEnergyModesThermo operator* <EquationOfState>
        (
            const scalar,
            const decoupledEnergyModesThermo&
        );

        friend decoupledEnergyModesThermo operator== <EquationOfState>
        (
            const decoupledEnergyModesThermo&,
            const decoupledEnergyModesThermo&
        );


    // Ostream Operator

        friend Ostream& operator<< <EquationOfState>
        (
            Ostream&,
            const decoupledEnergyModesThermo&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "decoupledEnergyModesThermoI.H"

#ifdef NoRepository
#   include "decoupledEnergyModesThermo.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
