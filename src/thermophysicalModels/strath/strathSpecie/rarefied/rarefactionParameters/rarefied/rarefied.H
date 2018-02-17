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
    Foam::rarefied

Description
    rarefied model.
  
SourceFiles
    rarefied.C

\*---------------------------------------------------------------------------*/

#ifndef rarefied_H
#define rarefied_H

#include "rarefactionParameter.H"

#include "fundamentalConstants.H"
#include "physicoChemicalConstants.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
  
/*---------------------------------------------------------------------------*\
                             Class rarefied Declaration
\*---------------------------------------------------------------------------*/

template<class ThermoType>
class rarefied
:
    public rarefactionParameter
{  

protected:
  
    // Protected data
    
        //- Thermodynamic data of the species
        const PtrList<ThermoType>& speciesThermo_;
        
        // Mole-fraction below which the contribution of the species is ignored
        const scalar miniXs_;


    // Protected member functions
    
        //- Update the coefficients of the species equations
        virtual void updateCoefficients();
    

public:

    //- Runtime type information
    TypeName("rarefied");

    
    // Constructors

        //- Construct from components
        rarefied
        (
            rho2ReactionThermo& thermo,
            const compressibleTurbulenceModel& turbulence
        );


      //- Destructor
      ~rarefied()
      {}


    // Member Functions

        //- Shear viscosity of a specie given by index.
        inline virtual scalar mu(label i, scalar p, scalar Tt) const
        {
            return speciesThermo_[i].mu(p, Tt);
        }
        
        //- Molecular weight [g/mol] of a specie given by index.
        inline scalar W(label i) const
        {
            return speciesThermo_[i].W();
        }
        
        //- Molecular weight ratio for a pair of species given by index.
        inline scalar ratioW(label i, label j) const
        {
            return W(i)/W(j);
        }
        
        //- Temperature exponent of viscosity for a specie given by index.
        inline scalar omega(label i) const
        {
            return speciesThermo_[i].omega();
        }
        
        //- Temperature-viscosity coefficient of a pair of species (i, j) 
        //  given by index.
        inline scalar omegaPair(label i, label j) const
        {
            if(mfpModelName_ == "hardSphere" or mfpModelName_ == "maxwellian")
            {
                return 0.5;
            }
            else
            {
                return 0.5*(omega(i)+omega(j));
            }
        }
        
        //- Diameter of a specie given by index.
        inline scalar diameter(label i) const
        {
            return speciesThermo_[i].diameter();
        }
        
        //- Diameter of a pair of species (i, j) given by index.
        inline scalar diameterPair(label i, label j) const
        {
            return 0.5*(diameter(i)+diameter(j));
        }
        
        // Calculations and updates the different fields
        virtual void correct(const volVectorField& U);
        
        //- Write the user-defined rarefied quantities
        virtual void write();
                
        //- Read transportProperties dictionary
        virtual bool read();
        
        //- Inner quantity calculation (Bird's book, 1993, Eq. 4.77)
        inline scalar innerQuantity(label i, label j, scalar Tt) const
        {
            return constant::mathematical::pi*sqr(diameterPair(i,j))
                *pow(273.0/Tt, omegaPair(i,j)-0.5)*sqrt(1.0 + ratioW(i,j));
        }
  };

  
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

}	// End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

#ifdef NoRepository
#   include "rarefied.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
