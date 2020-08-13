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
    Foam::Euler2Implicit

Description
    An Euler implicit solver for chemistry

SourceFiles
    Euler2Implicit.C

\*---------------------------------------------------------------------------*/

#ifndef Euler2Implicit_H
#define Euler2Implicit_H

#include "chemistry2Solver.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class Euler2Implicit Declaration
\*---------------------------------------------------------------------------*/

template<class Chemistry2Model>
class Euler2Implicit
:
    public chemistry2Solver<Chemistry2Model>
{
    // Private data

        //- Coefficients dictionary
        dictionary coeffsDict_;

        // Model constants

            //- Chemistry timescale
            scalar cTauChem_;

            //- Equilibrium rate limiter flag (on/off)
            Switch eqRateLimiter_;

        // Solver data
        mutable scalarField cTp_;


public:

    //- Runtime type information
    TypeName("Euler2Implicit");


    // Constructors

        //- Construct from mesh
        Euler2Implicit(const fvMesh& mesh);


    //- Destructor
    virtual ~Euler2Implicit();


    // Member Functions

        void updateRRInReactionI
        (
            const label index,
            const scalar pr,
            const scalar pf,
            const scalar corr,
            const label lRef,
            const label rRef,
            const scalar p,
            simpleMatrix<scalar>& RR
        ) const;

        //- Update the concentrations and return the chemical time
        virtual void solve
        (
            scalarField& c,
            scalar& T,
            List<scalar>& spTv,
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
            List<scalar>& spTv,
            scalar& p,
            scalar& deltaT,
            scalar& subDeltaT
        ) const;
        
        //- Update the concentrations and return the chemical time
        virtual void solve
        (
            scalarField& c,
            scalarField& cfwdiir,
            labelList& iirIds,
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
            scalarField& cfwdiir,
            labelList& iirIds,
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
#   include "Euler2Implicit.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
