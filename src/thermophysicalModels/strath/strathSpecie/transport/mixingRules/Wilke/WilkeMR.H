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
    Foam::WilkeMR

Description
    WilkeMR model.

SourceFiles
    WilkeMR.C

\*---------------------------------------------------------------------------*/

#ifndef WilkeMR_H
#define WilkeMR_H

#include "mixingRule.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                             Class WilkeMR Declaration
\*---------------------------------------------------------------------------*/

template<class ThermoType>
class WilkeMR
:
    public mixingRule
{

protected:

    // Protected data

        //- Thermodynamic data of the species
        const PtrList<ThermoType>& speciesThermo_;
        
        // Weight coefficient in WilkeMR's mixing rule
        PtrList<volScalarField> phi_;

        // Molar-fraction below which the contribution of the species is ignored
        const scalar miniXs_;


    // Protected member functions

        //- Update the coefficients of the species equations
        virtual void updatePhi();


public:

    //- Runtime type information
    TypeName("WilkeMR");


    // Constructors

        //- Construct from components
        WilkeMR
        (
            rho2ReactionThermo& thermo,
            const compressibleTurbulenceModel& turbulence
        );


      //- Destructor
      ~WilkeMR()
      {}


    // Member Functions

        //- Molar weight for a specie given by index
        virtual inline scalar W(const label i) const
        {
            return speciesThermo_[i].W();
        }

        //- Viscosity for a specie given by index
        virtual inline scalar mu
        (
            const label i,
            const scalar p,
            const scalar T
        ) const
        {
            return speciesThermo_[i].mu(p, T);
        }

        //- Thermal conductivity [W/mK], for a specie given by index
        //  trans-rotational energy mode
        virtual inline scalar kappatr
        (
            const label i,
            const scalar p,
            const scalar T
        ) const
        {
            return speciesThermo_[i].kappatr(p, T);
        }

        //- Thermal conductivity [W/mK], for a specie given by index
        //  vibro-electronic energy mode
        virtual inline scalar kappave
        (
            const label i,
            const scalar p,
            const scalar T,
            const scalar Tve
        ) const
        {
            return speciesThermo_[i].kappave(p, T, Tve);
        }

        // Calculations and updates of the different fields
        virtual void correct();

        //- Write the user-defined WilkeMR quantities
        virtual void write();

        //- Read transportProperties dictionary
        virtual bool read();
  };


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

}	// End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

#ifdef NoRepository
#   include "WilkeMR.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
