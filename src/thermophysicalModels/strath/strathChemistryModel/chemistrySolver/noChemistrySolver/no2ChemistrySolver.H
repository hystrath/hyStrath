/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

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
    Foam::no2ChemistrySolver

Description
    Dummy chemistry solver for 'none' option

SourceFiles
    no2ChemistrySolver.H
    no2ChemistrySolver.C

\*---------------------------------------------------------------------------*/

#ifndef noChemistySolver_H
#define noChemistySolver_H

#include "chemistry2Solver.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class no2ChemistrySolver Declaration
\*---------------------------------------------------------------------------*/

template<class Chemistry2Model>
class no2ChemistrySolver
:
    public chemistry2Solver<Chemistry2Model>
{

public:

    //- Runtime type information
    TypeName("none");


    // Constructors

        //- Construct from components
        no2ChemistrySolver(const fvMesh& mesh);


    //- Destructor
    virtual ~no2ChemistrySolver();


    // Member Functions

        //- Update the concentrations and return the chemical time
        virtual void solve
        (
            scalarField& c,
            scalar& T,
            scalar& Tv, // NEW VINCENT
            List<scalar>& spTv, // NEW VINCENT
            scalar& p,
            scalar& deltaT,
            scalar& subDeltaT
        ) const;
        
        //- Update the concentrations (also forward apart) and return the chemical time
        virtual void solve
        (
            scalarField& c,
            scalarField& cfwd,
            scalar& rriirN, // NEW VINCENT
            scalar& rriirO, // NEW VINCENT
            scalar& T,
            scalar& Tv, // NEW VINCENT
            List<scalar>& spTv, // NEW VINCENT
            scalar& p,
            scalar& deltaT,
            scalar& subDeltaT
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "no2ChemistrySolver.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
