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
    fourierPolyLeastSquaresFit

Description

\*----------------------------------------------------------------------------*/

#include "fourierPolyLeastSquaresFit.H"
#include "graph.H"
#include "simpleMatrix.H"

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //



fourierPolyLeastSquaresFit::fourierPolyLeastSquaresFit
(
    const scalarField& x,
    const scalarField& y,         
    const label& degree,
    const scalar& length
)
:
    coeffs_(degree, 0.0)
{
    
    // check if x and y are of the same size
    if(x.size() != y.size())
    {
        FatalErrorIn("fourierPolyLeastSquaresFit::fourierPolyLeastSquaresFit()")
            << "Error in input lists. x and y have to be of the same size. "
            << nl << " size of x: " << x.size() 
            << ", size of y: " << y.size()
            << exit(FatalError);
    }    
    
    // check that n >= m + 1
    if(x.size() < degree )
    {
        FatalErrorIn("fourierPolyLeastSquaresFit::fourierPolyLeastSquaresFit()")
            << "It is always appropriate that n >= k"
            << ", where n (= "<< x.size() << ") is the size of the input list (x,y)"
            << ", and k (= " << degree <<") is the degree of fourier polynomial."
            << exit(FatalError);
    }
    
    // check if n is odd
    
    // test if Nmicro are odd
    if(degree % 2 == 0)
    {
        FatalErrorIn("fourierPolyLeastSquaresFit::fourierPolyLeastSquaresFit()")
            << "The fourier polynomial is even = " << degree
            << ". It has to be odd."
            << exit(FatalError);               
    }    
    
    label matrixSize = degree;
    
    simpleMatrix<scalar> luMatrix(matrixSize, 0.0, 0.0);
    
    label N = (degree-1)*0.5;
    
    label row = 0;    
    
    // first row
    {
        
        scalar s0 = 0.0;
        
        forAll(x, n)
        {
            s0 += 1.0;
        }
        
        luMatrix[row][0] = s0;
        
        for (label j=1; j<= N; j++)
        {
            forAll(x, n)
            {
                luMatrix[row][(2*j)-1] += Foam::sin(2.0*constant::mathematical::pi*j*x[n]/length);
                luMatrix[row][(2*j)+1-1] += Foam::cos(2.0*constant::mathematical::pi*j*x[n]/length);        
            }
        }

        scalar tk = 0.0;
        
        forAll(y, n)
        {
            tk += y[n];
        }
        
        luMatrix.source()[row] = tk;
    }
    
    // for all other rows 
    
    for (label jdash=1; jdash<= N; jdash++)
    {
        row++;
        
        scalar s0 = 0.0;
        
        forAll(x, n)
        {
            s0 += Foam::sin(2.0*constant::mathematical::pi*jdash*x[n]/length);
        }
        
        luMatrix[row][0] = s0;
        
        // for all columns / terms in each equation
        for (label j=1; j<= N; j++)
        {
            forAll(x, n)
            {
                luMatrix[row][(2*j)-1] += Foam::sin(2.0*constant::mathematical::pi*j*x[n]/length)
                                          *Foam::sin(2.0*constant::mathematical::pi*jdash*x[n]/length);
                                          
                luMatrix[row][(2*j)+1-1] += Foam::cos(2.0*constant::mathematical::pi*j*x[n]/length)
                                            *Foam::sin(2.0*constant::mathematical::pi*jdash*x[n]/length);
            }
        }
    
        scalar tk = 0.0;
        
        forAll(y, n)
        {
            tk += y[n]*Foam::sin(2.0*constant::mathematical::pi*jdash*x[n]/length);
        }
        
        luMatrix.source()[row] = tk;
    }

    for (label jdash=1; jdash<= N; jdash++)
    {
        row++;
        
        scalar s0 = 0.0;
        
        forAll(x, n)
        {
            s0 += Foam::cos(2.0*constant::mathematical::pi*jdash*x[n]/length);
        }
        
        luMatrix[row][0] = s0;
        
        // for all columns / terms in each equation
        for (label j=1; j<= N; j++)
        {
            forAll(x, n)
            {
                luMatrix[row][(2*j)-1] += Foam::sin(2.0*constant::mathematical::pi*j*x[n]/length)
                                            *Foam::cos(2.0*constant::mathematical::pi*jdash*x[n]/length);
                                            
                luMatrix[row][(2*j)+1-1] += Foam::cos(2.0*constant::mathematical::pi*j*x[n]/length)
                                          *Foam::cos(2.0*constant::mathematical::pi*jdash*x[n]/length);
            }
        }
    
        scalar tk = 0.0;
        
        forAll(y, n)
        {
            tk += y[n]*Foam::cos(2.0*constant::mathematical::pi*jdash*x[n]/length);
        }
        
        luMatrix.source()[row] = tk;
    }

//     Info << "matrix = " << luMatrix << endl;
    
    coeffs_ = luMatrix.LUsolve();    
}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

fourierPolyLeastSquaresFit::~fourierPolyLeastSquaresFit()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const scalarField& fourierPolyLeastSquaresFit::coeffs() const
{
    return coeffs_;
}

// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //




// * * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
