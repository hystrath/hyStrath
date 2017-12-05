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
    linearLeastSquaresFit

Description

\*----------------------------------------------------------------------------*/

#include "linearLeastSquaresFit.H"
#include "graph.H"

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //
// void linearLeastSquaresFit::setRadius()
// {
//     for(label i = 0; i < noOfBins_; i++)
//     {
//        radius_[i] = (0.5 + scalar(i)) * binWidth();
//     }
// }

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


// y = ax + b
linearLeastSquaresFit::linearLeastSquaresFit
(
    const scalarField& x,
    const scalarField& y,         
//     scalarField sig,
//     label mwt,
    scalar& a,
    scalar& b
//     scalar& siga,
//     scalar& sigb,
//     scalar& chi2,
//     scalar& q
)

{

   
    a = 0.0;
    label nData = x.size();
    if(nData > 0)
    { 
    scalar sx= 0.0;
    scalar sy= 0.0;
    scalar ss = 0.0;
    scalar sxoss = 0.0;
    scalar t =0.0;
    scalar st2 =0.0;
    
    for (label i=0; i < nData; i++)
    {
        sx += x[i];
        sy += y[i];
    }
    
    ss = nData;
    
    sxoss=sx/ss;
    
    for (label i=0; i < nData; i++)
    {
        t = x[i]-sxoss;
        st2 += t*t;
        a += t*y[i];
    }
    
    a /= st2;
    b=(sy-(sx*a))/ss;
//     siga=sqrt((1.0+sx*sx/(ss*st2))/ss);
//     sigb=sqrt(1.0/st2);    
//     chi2=0.0;
//     q=0.0;
    }
    else
    {
	Info << "WARNING: least-squares fit has no data" << endl;
    }
}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

linearLeastSquaresFit::~linearLeastSquaresFit()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //




// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

// void linearLeastSquaresFit::operator=(const linearLeastSquaresFit& rhs)
// {
//     // Check for assignment to self
//     if (this == &rhs)
//     {
//         FatalErrorIn("linearLeastSquaresFit::operator=(const linearLeastSquaresFit&)")
//             << "Attempted assignment to self"
//             << abort(FatalError);
//     }
// 
//     Map<label>::operator=(rhs);
// 
//     binWidth_ = rhs.binWidth();
// }


// * * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

// Ostream& operator<<(Ostream& os, const linearLeastSquaresFit& d)
// {
//     os  << d.binWidth_
//         << static_cast<const Map<label>&>(d);
// 
//     // Check state of Ostream
//     os.check
//     (
//         "Ostream& operator<<(Ostream&, "
//         "const linearLeastSquaresFit&)"
//     );
// 
//     return os;
// }

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
