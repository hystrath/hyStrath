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

#include "Euler2Implicit.H"
#include "addToRunTimeSelectionTable.H"
#include "simpleMatrix.H"

#include "rho2ReactionThermo.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Chemistry2Model>
Foam::Euler2Implicit<Chemistry2Model>::Euler2Implicit
(
    const fvMesh& mesh
)
:
    chemistry2Solver<Chemistry2Model>(mesh),
    coeffsDict_(this->subDict("Euler2ImplicitCoeffs")),
    cTauChem_(readScalar(coeffsDict_.lookup("cTauChem"))),
    eqRateLimiter_(coeffsDict_.lookup("equilibriumRateLimiter")),
    cTp_(this->nEqns())
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Chemistry2Model>
Foam::Euler2Implicit<Chemistry2Model>::~Euler2Implicit()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Chemistry2Model>
void Foam::Euler2Implicit<Chemistry2Model>::updateRRInReactionI
(
    const label index,
    const scalar pr,
    const scalar pf,
    const scalar corr,
    const label lRef,
    const label rRef,
    const scalar p,
    simpleMatrix<scalar>& RR
) const
{
    const Reaction2<typename Chemistry2Model::thermoType>& R =
        this->reactions_[index];

    const label lSize = int(R.lhs().size());
    for(int s = 0; s < lSize; ++s)
    {
        label si = R.lhs()[s].index;
        scalar sl = R.lhs()[s].stoichCoeff;
        RR[si][rRef] -= sl*pr*corr;
        RR[si][lRef] += sl*pf*corr;
    }

    const label rSize = int(R.rhs().size());
    for(int s = 0; s < rSize; ++s)
    {
        label si = R.rhs()[s].index;
        scalar sr = R.rhs()[s].stoichCoeff;
        RR[si][lRef] -= sr*pf*corr;
        RR[si][rRef] += sr*pr*corr;
    }
}


template<class Chemistry2Model>
void Foam::Euler2Implicit<Chemistry2Model>::solve
(
    scalarField& c,
    scalar& T,
    scalarList& spTv,
    scalar& p,
    scalar& deltaT,
    scalar& subDeltaT
) const
{
    const label nSpecie = this->nSpecie();
    simpleMatrix<scalar> RR(nSpecie, 0, 0);

    for (label i=0; i<nSpecie; i++)
    {
        c[i] = max(0.0, c[i]);
    }

    scalar cTot = sum(c);

    scalar deltaTEst = min(deltaT, subDeltaT);

    forAll(this->reactions(), i)
    {
        scalar pf, cf, pr, cr;
        label lRef, rRef;

        scalar omegai =
            this->omegaI(i, c, T, spTv, p, pf, cf, lRef, pr, cr, rRef);

        scalar corr = 1.0;
        if (eqRateLimiter_)
        {
            if (omegai < 0.0)
            {
                corr = 1.0/(1.0 + pr*deltaTEst);
            }
            else
            {
                corr = 1.0/(1.0 + pf*deltaTEst);
            }
        }

        updateRRInReactionI(i, pr, pf, corr, lRef, rRef, p, RR);
    }

    // Calculate the stable/accurate time-step
    scalar tMin = GREAT;

    for (label i=0; i<nSpecie; i++)
    {
        scalar d = 0;
        for (label j=0; j<nSpecie; j++)
        {
            d -= RR[i][j]*c[j];
        }

        if (d < -SMALL)
        {
            tMin = min(tMin, -(c[i] + SMALL)/d);
        }
        else
        {
            d = max(d, SMALL);
            scalar cm = max(cTot - c[i], 1.0e-5);
            tMin = min(tMin, cm/d);
        }
    }
    
    subDeltaT = cTauChem_*tMin;
    deltaT = min(deltaT, subDeltaT);

    // Add the diagonal and source contributions from the time-derivative
    for (label i=0; i<nSpecie; i++)
    {
        RR[i][i] += 1.0/deltaT;
        RR.source()[i] = c[i]/deltaT;
    }

    // Solve for the new composition
    c = RR.LUsolve();

    // Limit the composition
    for (label i=0; i<nSpecie; i++)
    {
        c[i] = max(0.0, c[i]);
    }
}


template<class Chemistry2Model>
void Foam::Euler2Implicit<Chemistry2Model>::solve
(
    scalarField& c,
    scalarField& cfwd,
    scalar& T,
    scalarList& spTv,
    scalar& p,
    scalar& deltaT,
    scalar& subDeltaT
) const
{
    const label nSpecie = this->nSpecie();
    simpleMatrix<scalar> RR(nSpecie, 0, 0);
    simpleMatrix<scalar> RRfwd(nSpecie, 0, 0);

    for (label i=0; i<nSpecie; i++)
    {
        c[i] = max(0.0, c[i]);
    }

    scalar cTot = sum(c);

    scalar deltaTEst = min(deltaT, subDeltaT);

    forAll(this->reactions(), i)
    {
        scalar pf, cf, pr, cr;
        label lRef, rRef;

        scalar omegai =
            this->omegaI(i, c, T, spTv, p, pf, cf, lRef, pr, cr, rRef);

        scalar corr = 1.0;
        if (eqRateLimiter_)
        {
            if (omegai < 0.0)
            {
                corr = 1.0/(1.0 + pr*deltaTEst);
            }
            else
            {
                corr = 1.0/(1.0 + pf*deltaTEst);
            }
        }

        updateRRInReactionI(i, pr, pf, corr, lRef, rRef, p, RR);
        updateRRInReactionI(i, 0, pf, corr, lRef, rRef, p, RRfwd);
    }

    // Calculate the stable/accurate time-step
    scalar tMin = GREAT;

    for (label i=0; i<nSpecie; i++)
    {
        scalar d = 0;
        for (label j=0; j<nSpecie; j++)
        {
            d -= RR[i][j]*c[j];
        }

        if (d < -SMALL)
        {
            tMin = min(tMin, -(c[i] + SMALL)/d);
        }
        else
        {
            d = max(d, SMALL);
            scalar cm = max(cTot - c[i], 1.0e-5);
            tMin = min(tMin, cm/d);
        }
    }

    subDeltaT = cTauChem_*tMin;
    deltaT = min(deltaT, subDeltaT);
    
    // Add the diagonal and source contributions from the time-derivative
    for (label i=0; i<nSpecie; i++)
    {
        RR[i][i] += 1.0/deltaT;
        RR.source()[i] = c[i]/deltaT;

        RRfwd[i][i] += 1.0/deltaT;
        RRfwd.source()[i] = cfwd[i]/deltaT;
    }

    // Solve for the new composition
    c = RR.LUsolve();
    cfwd = RRfwd.LUsolve();

    // Limit the composition
    for (label i=0; i<nSpecie; i++)
    {
        c[i] = max(0.0, c[i]);
        cfwd[i] = max(0.0, cfwd[i]);
    }
}
// ************************************************************************* //
