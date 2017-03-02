
/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012 OpenFOAM Foundation
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
    Foam::sensible2Enthalpy

Description
    Thermodynamics mapping class to expose the sensible enthalpy function
    as the standard enthalpy function h(T).

\*---------------------------------------------------------------------------*/

#ifndef sensible2Enthalpy_H
#define sensible2Enthalpy_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class sensible2Enthalpy Declaration
\*---------------------------------------------------------------------------*/

template<class Thermo>
class sensible2Enthalpy
{

public:

    // Constructors

        //- Construct
        sensible2Enthalpy()
        {}


    // Member Functions

        //- Return the instantiated type name
        static word typeName()
        {
            return "sensible2Enthalpy";
        }

        // Fundamental properties

            static word name()
            {
                return "h";
            }

            // Sensible enthalpy [J/kmol]
            scalar he
            (
                const Thermo& thermo,
                const scalar p,
                const scalar Tt,
                const scalar Tv = 0
            ) const
            {
                return thermo.es(p, Tt, Tv); //TODO hs
            }
            
            // NEW VINCENT ****************************************************
            // Sensible enthalpy, trans-rotational mode [J/kmol]
            scalar het
            (
                const Thermo& thermo,
                const scalar p,
                const scalar Tt
            ) const
            {
                return thermo.ets(p, Tt); //TODO hts
            }
            
            // Sensible enthalpy, vibrational mode [J/kg], mode vibMode
            scalar hev_mode
            (
                const Thermo& thermo,
                const label vibMode,
                const scalar p,
                const scalar Tv
            ) const
            {
                return thermo.hvs_mode(vibMode, p, Tv);
            }
            
            // Sensible enthalpy, vibrational mode [J/kg]
            scalar hev
            (
                const Thermo& thermo,
                const scalar p,
                const scalar Tv
            ) const
            {
                return thermo.hvs(p, Tv);
            }
            
            // Sensible enthalpy, electronic mode [J/kg]
            scalar heel
            (
                const Thermo& thermo,
                const scalar p,
                const scalar Tv
            ) const
            {
                return thermo.hels(p, Tv);
            }
            
            // Sensible enthalpy, vib-electronic mode [J/kg], mode vibMode
            scalar hevel_mode
            (
                const Thermo& thermo,
                const label vibMode,
                const scalar p,
                const scalar Tv
            ) const
            {
                return thermo.hvels_mode(vibMode, p, Tv);
            }
            
            // Sensible enthalpy, vib-electronic mode [J/kg]
            scalar hevel
            (
                const Thermo& thermo,
                const scalar p,
                const scalar Tv
            ) const
            {
                return thermo.hvels(p, Tv);
            }
            // END NEW VINCENT ************************************************

            // Heat capacity at constant pressure [J/(kmol K)]
            scalar cpv
            (
                const Thermo& thermo,
                const scalar p,
                const scalar Tt,
                const scalar Tv = 0
            ) const
            {
                return thermo.cv(p, Tt, Tv); //TODO cp
            }

            //- cp/cp []
            scalar cpBycpv
            (
                const Thermo& thermo,
                const scalar p,
                const scalar Tt,
                const scalar Tv = 0
            ) const
            {
                return 1;
            }

            // Sensible enthalpy [J/kg]
            scalar HE
            (
                const Thermo& thermo,
                const scalar p,
                const scalar Tt,
                const scalar Tv = 0
            ) const
            {
                return thermo.Es(p, Tt, Tv); //TODO Hs
            }
            
            // NEW VINCENT ****************************************************
            // Sensible enthalpy, trans-rotational mode [J/kg]
            scalar HEt
            (
                const Thermo& thermo,
                const scalar p,
                const scalar Tt
            ) const
            {
                return thermo.Ets(p, Tt); //TODO Hts
            }
            
            // Sensible enthalpy, vibrational mode [J/kg], mode vibMode
            scalar HEv_mode
            (
                const Thermo& thermo,
                const label vibMode,
                const scalar p,
                const scalar Tv
            ) const
            {
                return thermo.Hvs_mode(vibMode, p, Tv);
            }
            
            // Sensible enthalpy, vibrational mode [J/kg]
            scalar HEv
            (
                const Thermo& thermo,
                const scalar p,
                const scalar Tv
            ) const
            {
                return thermo.Hvs(p, Tv);
            }
            
            // Sensible enthalpy, electronic mode [J/kg]
            scalar HEel
            (
                const Thermo& thermo,
                const scalar p,
                const scalar Tv
            ) const
            {
                return thermo.Hels(p, Tv);
            }
            
            // Sensible enthalpy, vib-electronic mode [J/kg], mode vibMode
            scalar HEvel_mode
            (
                const Thermo& thermo,
                const label vibMode,
                const scalar p,
                const scalar Tv
            ) const
            {
                return thermo.Hvels_mode(vibMode, p, Tv);
            }
            
            // Sensible enthalpy, vib-electronic mode [J/kg]
            scalar HEvel
            (
                const Thermo& thermo,
                const scalar p,
                const scalar Tv
            ) const
            {
                return thermo.Hvels(p, Tv);
            }
            // END NEW VINCENT ************************************************

            // DELETION VINCENT ***********************************************
            //- Temperature from sensible enthalpy
            //  given an initial temperature T0
            scalar THE
            (
                const Thermo& thermo,
                const scalar h,
                const scalar p,
                const scalar T0
            ) const
            {
                return thermo.TEs(h, p, T0); // TODO if ever used: thermo.THs(h, p, T0);
            }
            // END DELETION VINCENT *******************************************
            
            // NEW VINCENT ****************************************************
            //- Translational-rotational temperature from sensible enthalpy
            //  given an initial temperature T0
            scalar TtHEt
            (
                const Thermo& thermo,
                const scalar h,
                const scalar p,
                const scalar T0
            ) const
            {
                return thermo.TtEts(h, p, T0); //TODO TtHts
            }
            
            //- Vibrational temperature from sensible enthalpy
            //  given an initial temperature T0, mode vibMode
            scalar TvHEv_mode
            (
                const Thermo& thermo,
                const label vibMode,
                const scalar h,
                const scalar p,
                const scalar T0
            ) const
            {
                return thermo.TvEvs_mode(vibMode, h, p, T0); //TODO
            }
            
            //- Vibrational temperature from sensible enthalpy
            //  given an initial temperature T0
            scalar TvHEv
            (
                const Thermo& thermo,
                const scalar h,
                const scalar p,
                const scalar T0
            ) const
            {
                return thermo.TvHvs(h, p, T0);
            }
            
            //- Vib-electronic temperature from sensible enthalpy
            //  given an initial temperature T0, mode vibMode
            scalar TvelHEvel_mode
            (
                const Thermo& thermo,
                const label vibMode,
                const scalar h,
                const scalar p,
                const scalar T0
            ) const
            {
                return thermo.TvelEvels_mode(vibMode, h, p, T0); //TODO
            }
            
            //- Vib-electronic temperature from sensible enthalpy
            //  given an initial temperature T0
            scalar TvelHEvel
            (
                const Thermo& thermo,
                const scalar h,
                const scalar p,
                const scalar T0
            ) const
            {
                return thermo.TvelHvels(h, p, T0);
            }
            // END NEW VINCENT ************************************************
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
