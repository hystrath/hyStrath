/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Class
    Foam::molarMR

Description
    molarMR model.
  
SourceFiles
    molarMR.C

\*---------------------------------------------------------------------------*/

#ifndef molarMR_H
#define molarMR_H

#include "mixingRule.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
  
/*---------------------------------------------------------------------------*\
                             Class molarMR Declaration
\*---------------------------------------------------------------------------*/

template<class ThermoType>
class molarMR
:
    public mixingRule
{  

protected:
  
    // Protected data
    
        //- Thermodynamic data of the species
        const PtrList<ThermoType>& speciesThermo_;


    // Protected member functions
    
    

public:

    //- Runtime type information
    TypeName("molarMR");

    
    // Constructors

        //- Construct from components
        molarMR
        (
            rho2ReactionThermo& thermo,
            const compressibleTurbulenceModel& turbulence
        );


      //- Destructor
      ~molarMR()
      {}


    // Member Functions

        //- Molar weight for a specie given by index. Defined here
        //  becouse depends on thermoType template
        virtual inline scalar W(const label i) const
        {
            return speciesThermo_[i].W();
        }
        
        //- Viscosity for a specie given by index. Defined here
        //  becouse depends on thermoType template
        virtual inline scalar mu(const label i, const scalar p, const scalar Tt) const
        {
            return speciesThermo_[i].mu(p, Tt);
        }
        
        //- Thermal conductivity [W/mK], trans-rotational mode, 
        //  for a specie given by index.
        virtual inline scalar kappatr(const label i, const scalar p, const scalar Tt) const
        {
            return speciesThermo_[i].kappatr(p, Tt);
        }
        
        //- Thermal conductivity [W/mK], vib-electronic mode, 
        //  for a specie given by index.
        virtual inline scalar kappave(const label i, const scalar p, const scalar Tt, const scalar Tve) const
        {
            return speciesThermo_[i].kappave(p, Tt, Tve);
        }
        
        //- Thermal diffusivity [kg/ms], trans-rotational mode, 
        //  for a specie given by index.
        virtual inline scalar alphatr(const label i, const scalar p, const scalar Tt) const
        {
            return speciesThermo_[i].alphatr(p, Tt);
        }
        
        //- Thermal diffusivity [kg/ms], vib-electronic mode, 
        //  for a specie given by index.
        virtual inline scalar alphave(const label i, const scalar p, const scalar Tt, const scalar Tve) const
        {
            return speciesThermo_[i].alphave(p, Tt, Tve);
        }
        
        // Calculations and updates of the different fields
        virtual void correct();
        
        //- Write the user-defined molarMR quantities
        virtual void write();
                
        //- Read transportProperties dictionary
        virtual bool read();
  };

  
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

}	// End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

#ifdef NoRepository
#   include "molarMR.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
