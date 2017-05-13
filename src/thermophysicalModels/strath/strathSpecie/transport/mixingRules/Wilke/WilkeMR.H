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
    
        // Weight coefficient in WilkeMR's mixing rule
        PtrList<volScalarField> phi_;
        
        // Mole-fraction below which the contribution of the species is ignored
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
            const compressible::turbulenceModel2& turbulence
        );


      //- Destructor
      ~WilkeMR()
      {}


    // Member Functions

        //- Molar weight for a specie given by index. Defined here
        //  becouse depends on thermoType template
        inline virtual scalar W(const label i) const
        {
            return thermo_.composition().W(i);
        }
        
        //- Viscosity for a specie given by index. Defined here
        //  becouse depends on thermoType template
        inline virtual scalar mu(const label i, const scalar p, const scalar Tt) const
        {
            return thermo_.composition().mu(i, p, Tt);
        }
        
        //- Thermal conductivity, trans-rotational mode, 
        //  for a specie given by index.
        inline virtual scalar kappatr(const label i, const scalar p, const scalar Tt) const
        {
            return thermo_.composition().kappatr(i, p, Tt);
        }
        
        //- Thermal conductivity, vib-electronic mode, 
        //  for a specie given by index.
        inline virtual scalar kappave(const label i, const scalar p, const scalar Tt, const scalar Tve) const
        {
            return thermo_.composition().kappave(i, p, Tt, Tve);
        }
        
        //- Thermal diffusivity [kg/ms], trans-rotational mode, 
        //  for a specie given by index.
        inline virtual scalar alphatr(const label i, const scalar p, const scalar Tt) const
        {
            return thermo_.composition().alphatr(i, p, Tt);
        }
        
        //- Thermal diffusivity [kg/ms], vib-electronic mode, 
        //  for a specie given by index.
        inline virtual scalar alphave(const label i, const scalar p, const scalar Tt, const scalar Tve) const
        {
            return thermo_.composition().alphave(i, p, Tt, Tve);
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
