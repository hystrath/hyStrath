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
    Foam::AppletonBray

Description
    J. P. Appleton and K. N. C. Bray (1964) relaxation time model.
    "The Conservation equations for a non-equilibrium plasma"
    J. Fluid Mech., vol. 20, part 4, pp. 659-672

SourceFiles
    AppletonBray.C

\*---------------------------------------------------------------------------*/

#ifndef AppletonBray_H
#define AppletonBray_H

#include "relaxationTimeModelHE.H"
#include "constants.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class AppletonBray Declaration
\*---------------------------------------------------------------------------*/

template<class ThermoType>
class AppletonBray
:
    public relaxationTimeModelHE
{

protected:

    // Protected data

        //- Thermodynamic data of the species
        const PtrList<ThermoType>& speciesThermo_;

        //- Position of the electrons in the species table
        const label electronListPosition_;
        
        //- Collision cross-section for electron / heavy-particle interactions
        const scalar sigma_er_;

        //- Two constant factors appearing in the expression of the collision
        //  cross section for electron-ion interactions
        //  (cf. Scalabrin's PhD, 2007, eq. 2.67)
        scalar sigma_eIon_factor1_, sigma_eIon_factor2_;
        
        //- Constant factor appearing in the expression of the energy transfer
        //  rate between heavy-particles and electrons
        scalar Qhe_factor_;


    // Protected member functions


public:

    //- Runtime type information
    TypeName("AppletonBray");


    // Constructors

        //- Construct from components
        AppletonBray
        (
            rho2ReactionThermo& thermo,
            const compressibleTurbulenceModel& turbulence
        );


      //- Destructor
      ~AppletonBray()
      {}


    // Member Functions

        //- Molar weight for a specie given by index. Defined here
        //  becouse depend on thermoType template
        virtual scalar W(label i) const
        {
            return speciesThermo_[i].W()*1.0e-3;
        }

        // Correct the h-e relaxation times and h-e energy transfer
        virtual void correct();

        //- Read transportProperties dictionary
        virtual bool read();
  };


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

}	// End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

#ifdef NoRepository
#   include "AppletonBray.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
