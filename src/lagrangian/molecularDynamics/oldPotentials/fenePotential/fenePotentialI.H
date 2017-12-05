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
    distribution

Description

\*---------------------------------------------------------------------------*/

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Destructors * * * * * * * * * * * * * * * //



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline scalar fenePotential::fFene(const scalar& r)
{
//     scalar ra = r - dr_;
//     scalar rb = r + dr_;
// 
//     scalar F = -numericalDifferentiation
//     (
//         ra,
//         r,
//         rb,
//         Ufene(ra),
//         Ufene(r),
//         Ufene(rb)
//     );

    scalar F = 0;
    
    if(r < r0in_)
    {
        F = -k_*r/(1.0-((r/r0_)*(r/r0_))); 
    }
    else
    {
        F = -k_*r0in_/(1.0-((r0in_/r0_)*(r0in_/r0_)));
    }
    
    return F;
}
    
inline scalar fenePotential::force (const scalar& r)
{    
//     label k1 = label(r/dr_);
//     
//     label k2 = k1+1;
//     
//     scalar f = 0.0;
//     
//     if(k2 < U_.size())
//     {
//         f = ( (f_[k2]-f_[k1])*(r- (k1*dr_))/(k2*dr_ - k1*dr_) ) + f_[k1];
//     }
    scalar f = fFene(r); 
    
    return f;
}
    
inline scalar fenePotential::potentialEnergy (const scalar& r)
{    
//     label k1 = label(r/dr_);
//     
//     label k2 = k1+1;
//     
//     scalar U = 0.0;
//     
//     if(k2 < U_.size())
//     {
//         U = ( (U_[k2]-U_[k1])*(r- (k1*dr_))/(k2*dr_ - k1*dr_) ) + U_[k1];
//     }
    
    scalar U = Ufene(r);
    
    return U;
}

    
inline scalar fenePotential::Ufene(const scalar& r)
{
    scalar U = 0;
    
    if(r < r0in_)
    {
        U = -0.5*k_*r0_*r0_*Foam::log(1- ((r/r0_)*(r/r0_)));
    }
    else
    {
       U = -0.5*k_*r0_*r0_*Foam::log(1- ((r0in_/r0_)*(r0in_/r0_)));
    }    
    
    return U;
}

// inline scalar fenePotential::numericalDifferentiation
// (
//     const scalar& ra,
//     const scalar& r,
//     const scalar& rb,
//     const scalar& Ua,
//     const scalar& U,
//     const scalar& Ub
// )
// {
// 
//     scalar denominator = (ra - r)*(ra - rb)*(r - rb);
// 
//     scalar a1 =
//     (
//         rb*rb*(Ua - U) + ra*ra*(U - Ub) + r*r*(Ub - Ua)
//     ) / denominator;
// 
//     scalar a2 =
//     (
//         rb*(U - Ua) + r*(Ua - Ub) + ra*(Ub - U)
//     ) / denominator;
// 
//     return a1 + 2.0*a2*r;
// }
    
inline const scalar& fenePotential::k() const
{
    return k_;
}
    
inline const scalar& fenePotential::rCut() const
{
    return rCut_;
}

inline const scalar& fenePotential::rCutSqr() const
{
    return rCutSqr_;
}
inline const scalar& fenePotential::r0() const
{
    return r0_;
}
// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
