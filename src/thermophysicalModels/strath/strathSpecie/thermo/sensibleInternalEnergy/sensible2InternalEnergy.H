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
    Foam::sensible2InternalEnergy

Description
    Thermodynamics mapping class to expose the sensible internal energy function
    as the standard internal energy function e(T, Tv).

\*---------------------------------------------------------------------------*/

#ifndef sensible2InternalEnergy_H
#define sensible2InternalEnergy_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class sensible2InternalEnergy Declaration
\*---------------------------------------------------------------------------*/

template<class Thermo>
class sensible2InternalEnergy
{

public:

    // Constructors

        //- Construct
        sensible2InternalEnergy()
        {}


    // Member Functions

        //- Return the instantiated type name
        static word typeName()
        {
            return "sensible2InternalEnergy";
        }

        // Fundamental properties

            static word name()
            {
                return "e";
            }

            //- Sensible internal energy [J/kmol]
            scalar he
            (
                const Thermo& thermo,
                const scalar p,
                const scalar T,
                const scalar Tv = 0
            ) const
            {
                return thermo.es(p, T, Tv);
            }

            // Sensible internal energy, trans-rotational energy mode [J/kmol]
            scalar het
            (
                const Thermo& thermo,
                const scalar p,
                const scalar T
            ) const
            {
                return thermo.ets(p, T);
            }

            // Sensible internal energy, vibrational energy mode [J/kg]
            // for a given mode vibMode
            scalar hev_mode
            (
                const Thermo& thermo,
                const label vibMode,
                const scalar p,
                const scalar Tv
            ) const
            {
                return thermo.evs_mode(vibMode, p, Tv);
            }

            // Sensible internal energy, vibrational energy mode [J/kg]
            scalar hev
            (
                const Thermo& thermo,
                const scalar p,
                const scalar Tv
            ) const
            {
                return thermo.evs(p, Tv);
            }

            // Sensible internal energy, electronic energy mode [J/kg]
            scalar heel
            (
                const Thermo& thermo,
                const scalar p,
                const scalar Tv
            ) const
            {
                return thermo.eels(p, Tv);
            }

            // Sensible internal energy, vibro-electronic energy mode [J/kg]
            // for a given mode vibMode
            scalar hevel_mode // TODO should not exist
            (
                const Thermo& thermo,
                const label vibMode,
                const scalar p,
                const scalar Tv
            ) const
            {
                return thermo.evels_mode(vibMode, p, Tv);
            }

            // Sensible internal energy, vibro-electronic energy mode [J/kg]
            scalar hevel
            (
                const Thermo& thermo,
                const scalar p,
                const scalar Tv
            ) const
            {
                return thermo.evels(p, Tv);
            }

            //- Heat capacity at constant volume [J/(kmol K)]
            scalar cpv
            (
                const Thermo& thermo,
                const scalar p,
                const scalar T,
                const scalar Tv = 0
            ) const
            {
                return thermo.cv(p, T, Tv);
            }

            //- cp/cv []
            scalar cpBycpv
            (
                const Thermo& thermo,
                const scalar p,
                const scalar T,
                const scalar Tv = 0
            ) const
            {
                return thermo.gamma(p, T, Tv);
            }

            //- Sensible internal energy [J/kg]
            scalar HE
            (
                const Thermo& thermo,
                const scalar p,
                const scalar T,
                const scalar Tv = 0
            ) const
            {
                return thermo.Es(p, T, Tv);
            }

            // Sensible internal energy, trans-rotational energy mode [J/kg]
            scalar HEt
            (
                const Thermo& thermo,
                const scalar p,
                const scalar T
            ) const
            {
                return thermo.Ets(p, T);
            }

            // Sensible internal energy, vibrational energy mode [J/kg]
            // for a given mode vibMode
            scalar HEv_mode
            (
                const Thermo& thermo,
                const label vibMode,
                const scalar p,
                const scalar Tv
            ) const
            {
                return thermo.Evs_mode(vibMode, p, Tv);
            }

            // Sensible internal energy, vibrational energy mode [J/kg]
            scalar HEv
            (
                const Thermo& thermo,
                const scalar p,
                const scalar Tv
            ) const
            {
                return thermo.Evs(p, Tv);
            }

            // Sensible internal energy, electronic energy mode [J/kg]
            scalar HEel
            (
                const Thermo& thermo,
                const scalar p,
                const scalar Tv
            ) const
            {
                return thermo.Eels(p, Tv);
            }

            // Sensible internal energy, vibro-electronic energy mode [J/kg]
            // for a given mode vibMode
            // TODO should not exist
            scalar HEvel_mode
            (
                const Thermo& thermo,
                const label vibMode,
                const scalar p,
                const scalar Tv
            ) const
            {
                return thermo.Evels_mode(vibMode, p, Tv);
            }

            // Sensible internal energy, vibro-electronic energy mode [J/kg]
            scalar HEvel
            (
                const Thermo& thermo,
                const scalar p,
                const scalar Tv
            ) const
            {
                return thermo.Evels(p, Tv);
            }

            //- Temperature from sensible internal
            //  energy given an initial temperature T0
            scalar THE
            (
                const Thermo& thermo,
                const scalar e,
                const scalar p,
                const scalar T0
            ) const
            {
                return thermo.TEs(e, p, T0);
            }

            //- Trans-rotational temperature from sensible internal
            //  energy given an initial temperature T0
            scalar TtHEt
            (
                const Thermo& thermo,
                const scalar e,
                const scalar p,
                const scalar T0
            ) const
            {
                return thermo.TtEts(e, p, T0);
            }

            //- Vibrational temperature from sensible internal
            //  energy given an initial temperature T0, mode vibMode
            scalar TvHEv
            (
                const Thermo& thermo,
                const label vibMode,
                const scalar e,
                const scalar p,
                const scalar T0
            ) const
            {
                return thermo.TvEvs_mode(vibMode, e, p, T0);
            }

            //- Vibrational temperature from sensible internal
            //  energy given an initial temperature T0
            scalar TvHEv
            (
                const Thermo& thermo,
                const scalar e,
                const scalar p,
                const scalar T0
            ) const
            {
                return thermo.TvEvs(e, p, T0);
            }

            //- Vibro-electronic temperature from sensible internal
            //  energy given an initial temperature T0, for a given mode vibMode
            scalar TvelHEvel_mode
            (
                const Thermo& thermo,
                const label vibMode,
                const scalar e,
                const scalar p,
                const scalar T0
            ) const
            {
                return thermo.TvelEvels_mode(vibMode, e, p, T0);
            }

            //- vibro-electronic temperature from sensible internal
            //  energy given an initial temperature T0
            scalar TvelHEvel
            (
                const Thermo& thermo,
                const scalar e,
                const scalar p,
                const scalar T0
            ) const
            {
                return thermo.TvelEvels(e, p, T0);
            }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
