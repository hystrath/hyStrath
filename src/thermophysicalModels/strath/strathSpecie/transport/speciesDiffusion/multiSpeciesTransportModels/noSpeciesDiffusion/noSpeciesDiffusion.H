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
    Foam::noSpeciesDiffusion

Description
    noSpeciesDiffusion multispecies transport model.

SourceFiles
    noSpeciesDiffusion.C

\*---------------------------------------------------------------------------*/

#ifndef noSpeciesDiffusion_H
#define noSpeciesDiffusion_H

#include "multiSpeciesTransportModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class noSpeciesDiffusion Declaration
\*---------------------------------------------------------------------------*/

template<class ThermoType>
class noSpeciesDiffusion
:
    public multiSpeciesTransportModel
{

protected:

    // Protected data

        //- Thermodynamic data of the species
        const PtrList<ThermoType>& speciesThermo_;


    // Protected member functions

        //- Update the coefficients of the species equations
        virtual void updateCoefficients();


public:

    //- Runtime type information
    TypeName("noSpeciesDiffusion");


    // Constructors

        //- Construct from components
        noSpeciesDiffusion
        (
            rho2ReactionThermo& thermo,
            const compressible::turbulenceModel& turbulence
        );


      //- Destructor
      ~noSpeciesDiffusion()
      {}


    // Member Functions

        //- Molar weight [kg/mol] for a specie given by index.
        virtual scalar W(label i) const
        {
            return speciesThermo_[i].W()*1.0e-3;
        }

        //- Absolute enthalpy [J/kg] for a specie given by index and temperature
        virtual scalar ha(label i, scalar p, scalar Tt, scalar Tv) const
        {
            return speciesThermo_[i].Ha(p, Tt, Tv);
        }
        
        //- Sensible enthalpy [J/kg] for a specie given by index and temperature
        virtual scalar hs(label i, scalar p, scalar Tt, scalar Tv) const
        {
            return speciesThermo_[i].Hs(p, Tt, Tv);
        }

        //- Update the species diffusion coefficients
        virtual void correct();

        //- Read transportProperties dictionary
        virtual bool read();
  };


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

}	// End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

#ifdef NoRepository
#   include "noSpeciesDiffusion.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
