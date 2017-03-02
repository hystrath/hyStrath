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
    Foam::decoupledEnergyModesThermo

Description
    uEM tables based thermodynamics package templated
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

        static const int nCoeffs_ = 7;
        typedef FixedList<scalar, nCoeffs_> cvCoeffsArray;
        
        typedef DynamicList<scalar> dynamicArray; // NEW VINCENT 11/02/2016


private:

    // Private data

        // Temperature limits of applicability of functions
        scalar Tlow_, Thigh_;

        cvCoeffsArray decoupledCvCoeffs_;
        dynamicArray vibrationalList_, electronicList_;


    // Private Member Functions

        //- Check that input data is valid
        void checkInputData() const;

        //- Return the cv coefficients
        inline const cvCoeffsArray& cvCoeffs(const scalar T) const;
        
        //- Return the vibrational coefficients in the format degeneracy - temperature
        inline const dynamicArray& vibrationalCoeffs(const scalar T) const;
        
        //- Return the electronic coefficients in the format degeneracy - temperature
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
        inline decoupledEnergyModesThermo(const word&, const decoupledEnergyModesThermo&);


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
            
            // NEW VINCENT ****************************************************
            //- Heat capacity at constant volume for the 
            //  translational mode [J/(kmol K)]
            inline scalar cv_trans(const scalar p, const scalar Tt) const;
            
            //- Heat capacity at constant volume for the 
            //  rotational mode [J/(kmol K)]
            inline scalar cv_r(const scalar p, const scalar Tt) const;
            
            //- Heat capacity at constant volume for the 
            //  trans-rotational mode [J/(kmol K)]
            inline scalar cv_t(const scalar p, const scalar Tt) const;
            
            //- Heat capacity at constant volume for the 
            //  vibrational energy mode, for the mode vibMode [J/(kmol K)]
            inline scalar cv_v_mode(const label vibMode, const scalar p, const scalar Tv) const;
            
            //- Heat capacity at constant volume for the 
            //  vibrational mode [J/(kmol K)]
            inline scalar cv_v(const scalar p, const scalar Tv) const;
            
            //- Heat capacity at constant volume for the 
            //  electronic mode [J/(kmol K)]
            inline scalar cv_el(const scalar p, const scalar Tv) const;
            
            //- Heat capacity at constant volume for the 
            //  electron mode [J/(kmol K)]
            inline scalar cv_e(const scalar p, const scalar Tv) const;
            // END NEW VINCENT ************************************************

            //- Absolute enthalpy [J/kmol]
            inline scalar ha(const scalar p, const scalar Tt, const scalar Tv = 0) const;

            // DELETION VINCENT ***********************************************
            //- Sensible enthalpy [J/kmol]
            //inline scalar hs(const scalar p, const scalar Tt, const scalar Tv = 0) const;
            // END DELETION VINCENT *******************************************
            
            // NEW VINCENT ****************************************************
            //- Trans-rotational sensible internal energy [J/kmol]
            inline scalar ets(const scalar p, const scalar Tt) const;
            
            //- Rotational sensible internal energy [J/kmol]
            inline scalar ers(const scalar p, const scalar Tt) const;
            
            //- Vibrational sensible internal energy for vib. mode vibMode [J/kmol]
            inline scalar evs_mode(const label vibMode, const scalar p, const scalar Tv) const;
            
            //- Vibrational sensible internal energy [J/kmol]
            inline scalar evs(const scalar p, const scalar Tv) const;
            
            //- Electronic sensible internal energy [J/kmol]
            inline scalar eels(const scalar p, const scalar Tv) const;
            
             //- Electron sensible internal energy [J/kmol]
            inline scalar ees(const scalar p, const scalar Tv) const;
            
            //- Degrees of freedom associated with the rotational energy mode [-]
            inline scalar zetar(const scalar p, const scalar Tt, const scalar Tv) const;
            
            //- Degrees of freedom associated with the vibrational energy mode, mode vibMode [-]
            inline scalar zetav_mode(const label vibMode, const scalar p, const scalar Tv) const;
            
            //- Degrees of freedom associated with the vibrational energy mode [-]
            inline scalar zetav(const scalar p, const scalar Tv) const;
            
            //- Degrees of freedom associated with the electronic energy mode [-]
            inline scalar zetael(const scalar p, const scalar Tv) const;
            
            //- Degrees of freedom associated with the electron energy mode [-]
            inline scalar zetae(const scalar p, const scalar Tv) const;
            // END NEW VINCENT ************************************************

            //- Chemical enthalpy [J/kmol]
            inline scalar hc() const;

            //- Entropy [J/(kmol K)]
            inline scalar s(const scalar p, const scalar Tt, const scalar Tv = 0) const;


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
