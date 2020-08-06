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
    Foam::ode2

Description
    An ODE solver for chemistry

SourceFiles
    ode2.C

\*---------------------------------------------------------------------------*/

#ifndef ode2_H
#define ode2_H

#include "chemistry2Solver.H"
#include "ODESolver.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                            Class ode2 Declaration
\*---------------------------------------------------------------------------*/

template<class Chemistry2Model>
class ode2
:
    public chemistry2Solver<Chemistry2Model>
{
    // Private data

        dictionary coeffsDict_;
        autoPtr<ODESolver> ode2Solver_;

        // Solver data
        mutable scalarField cTp_;


public:

    //- Runtime type information
    TypeName("ode2");


    // Constructors

        //- Construct from mesh
        ode2(const fvMesh& mesh);


    //- Destructor
    virtual ~ode2();


    // Member Functions

        //- Update the concentrations and return the chemical time
        virtual void solve
        (
            scalarField& c,
            scalar& T,
            scalarList& spTv,
            scalar& p,
            scalar& deltaT,
            scalar& subDeltaT
        ) const;

        //- Update the concentrations and return the chemical time
        virtual void solve
        (
            scalarField& c,
            scalarField& cfwd,
            scalar& T,
            scalarList& spTv,
            scalar& p,
            scalar& deltaT,
            scalar& subDeltaT
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "ode2.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
