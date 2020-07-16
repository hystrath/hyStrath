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
    Foam::pairPotentials::FENE

Description


SourceFiles
    FENE.C

\*---------------------------------------------------------------------------*/

#ifndef FENE_H
#define FENE_H

#include "pairPotential.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace pairPotentials
{

/*---------------------------------------------------------------------------*\
                         Class FENE Declaration
\*---------------------------------------------------------------------------*/

class FENE
:
    public pairPotential
{
    // Private data

        dictionary FENECoeffs_;

        scalar k_;
        scalar rO_;


public:

    //- Runtime type information
    TypeName("FENE");


    // Constructors

        //- Construct from components
        FENE
        (
            const word& name,
            const reducedUnits& rU,
            const dictionary& pairPotentialProperties
        );


    // Destructor

        ~FENE()
        {}


    // Member Functions

        scalar unscaledEnergy(const scalar r) const;

//         const scalar& sigma() const;
//         const scalar& epsilon() const;

        const dictionary& dict() const;

        //- Read dictionary
        bool read
        (
            const dictionary& pairPotentialProperties,
            const reducedUnits& rU
        );


};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace pairPotentials
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
