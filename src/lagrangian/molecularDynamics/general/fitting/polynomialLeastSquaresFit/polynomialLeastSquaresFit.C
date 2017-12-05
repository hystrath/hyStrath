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
    polynomialLeastSquaresFit

Description

\*----------------------------------------------------------------------------*/

#include "polynomialLeastSquaresFit.H"
#include "graph.H"
#include "simpleMatrix.H"

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //



polynomialLeastSquaresFit::polynomialLeastSquaresFit
(
    const scalarField& x,
    const scalarField& y,         
    const label& degree
)
:
    coeffs_(degree+1, 0.0)
{
    
    // check if x and y are of the same size
    if(x.size() != y.size())
    {
        FatalErrorIn("polynomialLeastSquaresFit::polynomialLeastSquaresFit()")
            << "Error in input lists. x and y have to be of the same size. "
            << nl << " size of x: " << x.size() 
            << ", size of y: " << y.size()
            << exit(FatalError);
    }    
    
    // check that n >= m + 1
    if(x.size() < degree+1 )
    {
        FatalErrorIn("polynomialLeastSquaresFit::polynomialLeastSquaresFit()")
            << "It is always appropriate that n >= m+1"
            << ", where n (= "<< x.size() << ") is the size of the input list (x,y)"
            << ", and m (= " << degree <<") is the degree of polynomial."
            << exit(FatalError);
    }     
    
    label matrixSize = degree + 1;
    
    simpleMatrix<scalar> luMatrix(matrixSize, 0.0, 0.0);
    
    // for all rows / equations
    for (label i = 0; i <= degree; i++)
    {
        label row = i;

        // for all columns / terms in each equation
        for (label j = 0; j <= degree; j++)
        {
            label col = j;
            
            scalar sk = 0.0;
            
            forAll(x, n)
            {
                sk += Foam::pow(x[n], (j+i) );
            }
            
            luMatrix[row][col] = sk;
        }
        
        scalar tk = 0.0;        
        
        forAll(x, n)
        {
            tk += y[n]*Foam::pow(x[n], i);
        }        
        
        luMatrix.source()[row] = tk; 
    }
    
//     Info << "matrix = " << luMatrix << endl;
    
    coeffs_ = luMatrix.LUsolve();    
}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

polynomialLeastSquaresFit::~polynomialLeastSquaresFit()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const scalarField& polynomialLeastSquaresFit::coeffs() const
{
    return coeffs_;
}

// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //




// * * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
