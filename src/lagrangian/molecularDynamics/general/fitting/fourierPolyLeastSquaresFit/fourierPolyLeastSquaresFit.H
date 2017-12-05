/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2005 OpenCFD Ltd.
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


Description
    Obtained from:

    http://mathworld.wolfram.com/LeastSquaresFittingPolynomial.html
    
    Solving for the polynomial equation:
    
    y = a_0 + a_1*x + ... + a_m x^m
    
    where m = degree.
    
    return coeffs in this order a_0, a_1, ... , a_m 
    
    (note the list is of size m + 1)
    
    ensure that n >= m + 1, wher n = size of input lists.
    
    
SourceFiles
    fourierPolyLeastSquaresFit.C

\*---------------------------------------------------------------------------*/

#ifndef fourierPolyLeastSquaresFit_H
#define fourierPolyLeastSquaresFit_H

#include "scalarField.H"
#include "polyMesh.H"
#include "IFstream.H"
#include "OFstream.H"
#include "Time.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class fourierPolyLeastSquaresFit Declaration
\*---------------------------------------------------------------------------*/

class fourierPolyLeastSquaresFit
{
    // Private data
    
         scalarField coeffs_;

    // Private Member Functions


public:

    // Static data members

    // Constructors

        //- Construct null

        //- Construct from xs and ys
        fourierPolyLeastSquaresFit
        (
            const scalarField& x,
            const scalarField& y, 
            const label& degree,
            const scalar& length
        );


    // Selectors

    // Destructor

        ~fourierPolyLeastSquaresFit();

    // Member Functions

        const scalarField& coeffs() const;

    // Member Operators

    // Friend Functions

    // Friend Operators

    // IOstream Operators

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
