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

\*---------------------------------------------------------------------------*/

#include "ode2.H"
#include "chemistry2Model.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Chemistry2Model>
Foam::ode2<Chemistry2Model>::ode2
(
    const fvMesh& mesh
)
:
    chemistry2Solver<Chemistry2Model>(mesh),
    coeffsDict_(this->subDict("ode2Coeffs")),
    ode2Solver_(ODESolver::New(*this, coeffsDict_)),
    cTp_(this->nEqns())
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Chemistry2Model>
Foam::ode2<Chemistry2Model>::~ode2()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Chemistry2Model>
void Foam::ode2<Chemistry2Model>::solve
(
    scalarField& c,
    scalar& T,
    scalarList& spTv,
    scalar& p,
    scalar& deltaT,
    scalar& subDeltaT
) const
{
    label nSpecie = this->nSpecie();

    // Copy the concentration, T, Tv, and p to the total solve-vector
    for (register int i=0; i<nSpecie; i++)
    {
        cTp_[i] = c[i];
    }
    cTp_[nSpecie] = T;
    cTp_[nSpecie+1] = spTv[0]; // NEW VINCENT
    cTp_[nSpecie+2] = p; // MODIFIED VINCENT

    ode2Solver_->solve(0, deltaT, cTp_, subDeltaT);

    for (register int i=0; i<nSpecie; i++)
    {
        c[i] = max(0.0, cTp_[i]);
    }
    T = cTp_[nSpecie];
    spTv[0] = cTp_[nSpecie+1]; // NEW VINCENT TODO WRONG REVISE IF NEEDED
    p = cTp_[nSpecie+2]; // MODIFIED VINCENT
}


template<class Chemistry2Model>
void Foam::ode2<Chemistry2Model>::solve
(
    scalarField&,
    scalarField&,
    scalar&,
    scalarList&,
    scalar&,
    scalar&,
    scalar&
) const
{}


// ************************************************************************* //
