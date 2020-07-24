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
    Foam::binaryDiffusionModels::Stephani

Description
    Stephani's binary diffusion model

    (see Development of a Hybrid DSMC/CFD Method for Hypersonic Boundary Layer
     Flow over Discrete Surface Roughness by Kelly Ann Stephani (May 2012)
     THE UNIVERSITY OF TEXAS at AUSTIN, Appendix C)

SourceFiles
    Stephani.C

\*---------------------------------------------------------------------------*/

#ifndef Stephani_H
#define Stephani_H

#include "binaryDiffusionModel.H"
#include "dimensionedScalar.H"
#include "volFields.H"

#include "constants.H"
#include <math.h>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace binaryDiffusionModels
{

/*---------------------------------------------------------------------------*\
                           Class Stephani Declaration
\*---------------------------------------------------------------------------*/

class Stephani
:
    public binaryDiffusionModel
{
    // Private data

        //- Fundamental constants (more efficient)
        const scalar pi, kB, Runi;

        //- Reduced molecular weight and reduced mass
        const scalar MStar_, mStar_;

        //- Parameters for Stephani's model
        const scalar omegaref_, dref_, Tref_, constantOmegaNeutral1_;


    // Private member functions

        inline scalar constantExpressionInOmegaNeutral1() const
        {
            return pi*sqr(dref_)/2.0*pow(Tref_, omegaref_ - 0.5)
                *sqrt(kB/(2.0*pi*mStar_))*tgamma(3.5-omegaref_)
                /tgamma(2.5-omegaref_);
        }

        inline scalar omegaNeutral1(const scalar T) const
        {
            return constantOmegaNeutral1_*pow(T, 1.0-omegaref_);
        }


public:

    //- Runtime type information
    TypeName("Stephani");


    // Constructors

        //- construct from components
        Stephani
        (
            const word& name1,
            const word& name2,
            const dictionary& dictThermo,
            const dictionary& dictTransport,
            const volScalarField& p,
            const volScalarField& pe,
            const volScalarField& T
        );


    // Destructor

        ~Stephani()
        {}


    // Member Functions

        //- Binary diffusion coefficient [m^2/s]
        virtual tmp<volScalarField> D() const;

        //- Binary diffusion coefficient for patch [m^2/s]
        virtual tmp<scalarField> D
        (
            const scalarField& p,
            const scalarField& pe,
            const scalarField& T,
            const label patchi
        ) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace binaryDiffusionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
