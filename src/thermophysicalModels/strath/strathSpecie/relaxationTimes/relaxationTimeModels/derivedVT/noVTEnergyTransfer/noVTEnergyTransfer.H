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
    Foam::noVTEnergyTransfer

Description
    noVTEnergyTransfer disables V-T energy transfer.

SourceFiles
    noVTEnergyTransfer.C

\*---------------------------------------------------------------------------*/

#ifndef noVTEnergyTransfer_H
#define noVTEnergyTransfer_H

#include "relaxationTimeModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class noVTEnergyTransfer Declaration
\*---------------------------------------------------------------------------*/

template<class ThermoType>
class noVTEnergyTransfer
:
    public relaxationTimeModel
{

protected:

    // Protected data

        //- Thermodynamic data of the species
        const PtrList<ThermoType>& speciesThermo_;


    // Protected member functions



public:

    //- Runtime type information
    TypeName("noVTEnergyTransfer");


    // Constructors

        //- Construct from components
        noVTEnergyTransfer
        (
            rho2ReactionThermo& thermo,
            const compressibleTurbulenceModel& turbulence
        );


      //- Destructor
      ~noVTEnergyTransfer()
      {}


    // Member Functions

        //- Molar weight [kg/mol] for a specie given by index. Defined here
        //  becouse depend on thermoType template
        virtual scalar W(label i) const
        {
            return speciesThermo_[i].W()*1.0e-3;
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
#   include "noVTEnergyTransfer.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
