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
    Foam::chemistry2Solver

Description
    An abstract base class for solving chemistry

SourceFiles
    chemistry2Solver.C

\*---------------------------------------------------------------------------*/

#ifndef chemistry2Solver_H
#define chemistry2Solver_H

#include "chemistry2Model.H"
#include "IOdictionary.H"
#include "scalarField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class chemistrySolver Declaration
\*---------------------------------------------------------------------------*/

template<class Chemistry2Model>
class chemistry2Solver
:
    public Chemistry2Model
{

public:

    // Constructors

        //- Construct from components
        chemistry2Solver(const fvMesh& mesh);


    //- Destructor
    virtual ~chemistry2Solver();


    // Member Functions

        //- Update the concentrations and return the chemical time
        virtual void solve
        (
            scalarField &c,
            scalar& T,
            scalar& Tv, // NEW VINCENT
            List<scalar>& spTv, // NEW VINCENT
            scalar& p,
            scalar& deltaT,
            scalar& subDeltaT
        ) const = 0;
        
        //- Update the concentrations (also forward apart) and return the chemical time
        virtual void solve
        (
            scalarField& c,
            scalarField& cwfd, // NEW VINCENT
            scalar& rriirN, // NEW VINCENT
            scalar& rriirO, // NEW VINCENT
            scalar& T,
            scalar& Tv, // NEW VINCENT
            List<scalar>& spTv, // NEW VINCENT
            scalar& p,
            scalar& deltaT,
            scalar& subDeltaT
        ) const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "chemistry2Solver.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
