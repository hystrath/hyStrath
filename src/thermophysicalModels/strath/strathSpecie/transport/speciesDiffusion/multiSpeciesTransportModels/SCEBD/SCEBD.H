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
    Foam::SCEBD

Description
    SCEBD multispecies transport model.

    Reference paper:
    J.D. Ramshaw and C.H. Chang,
    "Friction-Weighted Self-Consistent Effective Binary Diffusion Approximation,"
    J. Non-Equilib. Thermodyn. 21, 223 (1996)

    http://pdxscholar.library.pdx.edu/cgi/viewcontent.cgi?article=1095&context=phy_fac

SourceFiles
    SCEBD.C

\*---------------------------------------------------------------------------*/

#ifndef SCEBD_H
#define SCEBD_H

#include "multiSpeciesTransportModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                             Class SCEBD Declaration
\*---------------------------------------------------------------------------*/

template<class ThermoType>
class SCEBD
:
    public multiSpeciesTransportModel
{

protected:

    // Protected data

        //- Thermodynamic data of the species
        const PtrList<ThermoType>& speciesThermo_;

        // Molar-fraction below which the contribution of the species is ignored
        const scalar miniXs_;


    // Protected member functions

        //- Update the coefficients of the species equations
        virtual void updateCoefficients();


public:

    //- Runtime type information
    TypeName("SCEBD");


    // Constructors

        //- Construct from components
        SCEBD
        (
            rho2ReactionThermo& thermo,
            const compressible::turbulenceModel& turbulence
        );


      //- Destructor
      ~SCEBD()
      {}


    // Member Functions

        //- Molar weight [kg/mol] for a specie given by index
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
#   include "SCEBD.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
