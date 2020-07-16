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
    Foam::LewisNumber

Description
    Multispecies Lewis number transport model.

SourceFiles
    LewisNumber.C

\*---------------------------------------------------------------------------*/

#ifndef LewisNumber_H
#define LewisNumber_H

#include "Fick.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                             Class LewisNumber Declaration
\*---------------------------------------------------------------------------*/

template<class ThermoType>
class LewisNumber
:
    public Fick<ThermoType>
{

protected:

    // Protected data

        //- Lewis number
        scalar Le_;

    // Protected member functions

        //- Update the coefficients of the species equations
        void updateCoefficients();


public:

    //- Runtime type information
    TypeName("LewisNumber");


    // Constructors

        //- Construct from components
        LewisNumber
        (
            rho2ReactionThermo& thermo,
            const compressible::turbulenceModel& turbulence
        );


    //- Destructor
    virtual ~LewisNumber()
    {}

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

}	// End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

#ifdef NoRepository
#   include "LewisNumber.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
