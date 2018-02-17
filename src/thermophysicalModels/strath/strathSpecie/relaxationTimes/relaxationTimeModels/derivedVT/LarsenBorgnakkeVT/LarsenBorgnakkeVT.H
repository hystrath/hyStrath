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
    Foam::LarsenBorgnakkeVT

Description
    LarsenBorgnakkeVT relaxation time model.
  
SourceFiles
    LarsenBorgnakkeVT.C

\*---------------------------------------------------------------------------*/

#ifndef LarsenBorgnakkeVT_H
#define LarsenBorgnakkeVT_H

#include "relaxationTimeModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
  
/*---------------------------------------------------------------------------*\
                             Class LarsenBorgnakkeVT Declaration
\*---------------------------------------------------------------------------*/

template<class ThermoType>
class LarsenBorgnakkeVT
:
    public relaxationTimeModel
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
    TypeName("LarsenBorgnakkeVT");

    
    // Constructors

        //- Construct from components
        LarsenBorgnakkeVT
        (
            rho2ReactionThermo& thermo,
            const compressibleTurbulenceModel& turbulence
        );


      //- Destructor
      ~LarsenBorgnakkeVT()
      {}


    // Member Functions

        //- Molar weight for a specie given by index. Defined here
        //  becouse depend on thermoType template
        scalar W(label i) const
        {
            return speciesThermo_[i].W();
        }
        
        //- Temperature exponent of viscosity for a specie given by index.
        //  Defined here becouse depend on thermoType template
        scalar omega(label i) const
        {
            return speciesThermo_[i].omega();
        }
        
        //- Reference diameter for a specie given by index.
        //  Defined here becouse depend on thermoType template
        scalar dref(label i) const
        {
            return speciesThermo_[i].diameter();
        }
        
        //- Reduced molecular mass for 2 species given by index.
        //  Defined here becouse depend on thermoType template
        scalar reducedMolecularMass(label i, label j) const
        {
            return W(i)*W(j)/(W(i)+W(j))/(6.02214129e26);
        }
        
        // Correct the V-T relaxation times and V-T energy transfer
        virtual void correct();
                
        //- Read transportProperties dictionary
        virtual bool read();
  };

  
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

}	// End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

#ifdef NoRepository
#   include "LarsenBorgnakkeVT.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
