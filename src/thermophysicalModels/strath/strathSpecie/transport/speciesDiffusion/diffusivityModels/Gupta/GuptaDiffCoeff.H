/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
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
    Foam::binaryDiffusivityModels::GuptaDiffCoeff

Description
    GuptaDiffCoeff binary diffusivity model (see Gupta, Roop N., Yos, Jerrold M. Thompson, Richard A. 1989).
    NASA Technical Memorandum 101528
    
SourceFiles
    GuptaDiffCoeff.C

\*---------------------------------------------------------------------------*/

#ifndef GuptaDiffCoeff_H
#define GuptaDiffCoeff_H

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
                           Class GuptaDiffCoeff Declaration
\*---------------------------------------------------------------------------*/

class GuptaDiffCoeff
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
        
        // Coeffs given in cm2*atm/sec
        inline scalar DijBar(const scalar T) const 
        {
            return 1.01325e1*exp(piOmega_[3])*pow(T, piOmega_[0]
                *sqr(log(T)) + piOmega_[1]*log(T) + piOmega_[2]);
        }    


public:

    //- Runtime type information
    TypeName("GuptaDiffCoeff");


    // Constructors

        //- construct from components
        GuptaDiffCoeff
        (
            const word& name1,
            const word& name2,
            const dictionary& dictThermo,
            const dictionary& dictTransport,
            const volScalarField& p,
            const volScalarField& T
        );


    // Destructor

        ~GuptaDiffCoeff()
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
