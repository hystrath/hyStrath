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
    Foam::LeeLandauTellereV

Description
    LeeLandauTellereV relaxation time model.

SourceFiles
    LeeLandauTellereV.C

\*---------------------------------------------------------------------------*/

#ifndef LeeLandauTellereV_H
#define LeeLandauTellereV_H

#include "relaxationTimeModeleV.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                             Class LeeLandauTellereV Declaration
\*---------------------------------------------------------------------------*/

template<class ThermoType>
class LeeLandauTellereV
:
    public relaxationTimeModeleV
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
    TypeName("LeeLandauTeller");


    // Constructors

        //- Construct from components
        LeeLandauTellereV
        (
            rho2ReactionThermo& thermo,
            const compressibleTurbulenceModel& turbulence
        );


      //- Destructor
      ~LeeLandauTellereV()
      {}


    // Member Functions

        //- Molar weight for a specie given by index. Defined here
        //  becouse depend on thermoType template
        scalar W(label i) const
        {
            return speciesThermo_[i].W();
        }

        // Correct the e-V relaxation times and e-V energy transfer
        virtual void correct();

        //- Read transportProperties dictionary
        virtual bool read();
  };


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

}	// End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

#ifdef NoRepository
#   include "LeeLandauTellereV.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
