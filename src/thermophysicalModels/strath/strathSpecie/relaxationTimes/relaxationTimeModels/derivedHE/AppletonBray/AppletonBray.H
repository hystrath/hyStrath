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
    Foam::AppletonBray

Description
    AppletonBray (1964) relaxation time model.
  
SourceFiles
    AppletonBray.C

\*---------------------------------------------------------------------------*/

#ifndef AppletonBray_H
#define AppletonBray_H

#include "relaxationTimeModelHE.H"

#include "fundamentalConstants.H"
#include "physicoChemicalConstants.H"

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
        label electronListPosition_;
        
        //- Fundamental and physico-chemical constants
        const scalar RR, NA, kB, ec, pi;


    // Protected member functions
    

public:

    //- Runtime type information
    TypeName("AppletonBray");

    
    // Constructors

        //- Construct from components
        AppletonBray
        (
            rho2ReactionThermo& thermo,
            const compressible::turbulenceModel2& turbulence
        );


      //- Destructor
      ~AppletonBray()
      {}


    // Member Functions

        //- Molar weight for a specie given by index. Defined here
        //  becouse depend on thermoType template
        virtual scalar W(label i) const
        {
            return speciesThermo_[i].W();
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
