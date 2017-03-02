/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
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
    Foam::binaryDiffusivityModels::Hirschfelder

Description
    Hirschfelder binary diffusivity model (see Sutton and Gnoffo 1998).
    
SourceFiles
    Hirschfelder.C

\*---------------------------------------------------------------------------*/

#ifndef Hirschfelder_H
#define Hirschfelder_H

#include "binaryDiffusivityModel.H"
#include "dimensionedScalar.H"
#include "volFields.H"

#include "constants.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace binaryDiffusivityModels
{

/*---------------------------------------------------------------------------*\
                           Class Hirschfelder Declaration
\*---------------------------------------------------------------------------*/

class Hirschfelder
:
    public binaryDiffusivityModel
{
    // Private data

        //- Coefficients of the model
        FixedList<scalar, 4> piOmega_;
        
        //- Molecular weights
        const scalar W1_, W2_;
        
        //- Fundamental constants (more efficient)
        const scalar pi, kB, Runi;

        
    // Private member functions
        
        //- Collision terms
        inline scalar collisionIntegralNeutral1(const scalar T) const 
        {
            if(T != 0.0)
            {
                return 8.0/3.0*sqrt(2.0*W1_*W2_/(pi*Runi*T*(W1_+W2_)))
                    *1.0e-20*piOmegaNeutral(T);
            }
            else
            {
                return 0.0;
            }
        }
        
        inline scalar piOmegaNeutral(const scalar T) const 
        {
            return piOmega_[3]*pow(T, piOmega_[0]
                *sqr(log(T)) + piOmega_[1]*log(T) + piOmega_[2]);
        }    


public:

    //- Runtime type information
    TypeName("Hirschfelder");


    // Constructors

        //- construct from components
        Hirschfelder
        (
            const word& name1,
            const word& name2,
            const dictionary& dictThermo,
            const dictionary& dictTransport,
            const volScalarField& p,
            const volScalarField& T
        );


    // Destructor

        ~Hirschfelder()
        {}


    // Member Functions

        //- Binary diffusion coefficient [m^2/s]
        virtual tmp<volScalarField> D() const;
         
        //- Binary diffusion coefficient for patch [m^2/s]
        virtual tmp<scalarField> D
        (
            const scalarField& p,
            const scalarField& T,
            const label patchi
        ) const;
         
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace binaryDiffusivityModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
