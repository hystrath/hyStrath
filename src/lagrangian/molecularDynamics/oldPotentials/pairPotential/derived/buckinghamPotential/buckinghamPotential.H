/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2008-2009 OpenCFD Ltd.
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
    Foam::pairPotentials::buckinghamPotential

Description


SourceFiles
    buckinghamPotential.C

\*---------------------------------------------------------------------------*/

#ifndef buckinghamPotential_H
#define buckinghamPotential_H

#include "pairPotential.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace pairPotentials
{

/*---------------------------------------------------------------------------*\
                         Class buckinghamPotential Declaration
\*---------------------------------------------------------------------------*/

class buckinghamPotential
:
    public pairPotential
{
    // Private data

        dictionary buckinghamPotentialCoeffs_;

        scalar rho_;
        scalar alpha_;
        scalar C_;        


public:

    //- Runtime type information
    TypeName("buckinghamPotential");


    // Constructors

        //- Construct from components
        buckinghamPotential
        (
            const word& name,
            const reducedUnits& rU,
            const dictionary& pairPotentialProperties
        );


    // Destructor

        ~buckinghamPotential()
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
