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

#include "simpleMatrix.H"

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Destructors * * * * * * * * * * * * * * * //



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
inline const scalar& brennerPotential::rCutMax() const
{
    return rIJ2_;
}


inline scalar brennerPotential::forceRepulsive(const scalar& r)
{
    scalar force = fRep(r);

    return force;
}

inline scalar brennerPotential::forceAttractive(const scalar& r)
{
    scalar force = fAtt(r);

    return force;
}

/*
inline scalar brennerPotential::force(const scalar& r, const scalar& Bij)
{
    scalar ra = r - dr_;
    scalar rb = r + dr_;

    scalar dFdrij = numericalDifferentiation
    (
        ra,
        r,
        rb,
        fIJ(ra),
        fIJ(r),
        fIJ(rb)
    );

    scalar dURepdrij = numericalDifferentiation
    (
        ra,
        r,
        rb,
        Urep(ra),
        Urep(r),
        Urep(rb)
    );

    scalar dUAttdrij = numericalDifferentiation
    (
        ra,
        r,
        rb,
        Uatt(ra),
        Uatt(r),
        Uatt(rb)
    );


    scalar F = -constantBrenner_ *
                (
                    dFdrij*Urep(r) +
                    dURepdrij*fIJ(r) -
                    Bij *
                    (
                        dFdrij*Uatt(r) + dUAttdrij*fIJ(r)
                    )
                );

    return F;
}
*/

// inline scalar brennerPotential::potentialEnergy(const scalar& r, const scalar& Bij)
// {
//     return fIJ(r)*constantBrenner_*(Urep(r) - Bij*Uatt(r));
// }



inline scalar brennerPotential::Urep(const scalar& r)
{
    return constantBrenner_*fIJ(r)*exp(constantRepulsive_*(r-rIJ0_));
}

inline scalar brennerPotential::Uatt(const scalar& r)
{
    return constantBrenner_*fIJ(r)*sIJ_*exp(constantAttractive_*(r-rIJ0_));
}



inline scalar brennerPotential::fRep(const scalar& r)
{
    scalar ra = r - dr_;
    scalar rb = r + dr_;

    scalar F = -numericalDifferentiation
    (
        ra,
        r,
        rb,
        Urep(ra),
        Urep(r),
        Urep(rb)
    );

    return F;
}

inline scalar brennerPotential::fAtt(const scalar& r)
{
    scalar ra = r - dr_;
    scalar rb = r + dr_;

    scalar F = -numericalDifferentiation
    (
        ra,
        r,
        rb,
        Uatt(ra),
        Uatt(r),
        Uatt(rb)
    );

    return F;
}



inline scalar brennerPotential::fIJ(const scalar& r)
{
    scalar fij = 1.0;

    if ( (r < rIJ2_) && (r >= rIJ1_) )
    {
        
       fij = 0.5*
             (
                1.0+
                cos
                (
                    constant::mathematical::pi*(r - rIJ1_)
                        /(rIJ2_-rIJ1_)
                )
             );
    }

    if(r >= rIJ2_)
    {
        fij = 0.0;
    }

    return fij;
}

inline scalar brennerPotential::F(const scalar& x)
{
    scalar f = 1.0;

    if ( (x < 3.0) && (x > 2.0) )
    {
       f = 0.5*(1.0+cos(constant::mathematical::pi*(x-2.0)));
    }

    if(x >= 3.0)
    {
        f = 0.0;
    }

    return f;
}

inline scalar brennerPotential::G(const scalar& cosThetaIJK)
{
    return a0_*(
                    1.0 + ((c0_*c0_)/(d0_*d0_))
                    -(
                        (c0_*c0_)/
                        (
                            (d0_*d0_) 
                            + (1.0 +cosThetaIJK)*(1.0+cosThetaIJK)
                        )
                    )
                );
}

inline scalar brennerPotential::numericalDifferentiation
(
    const scalar& ra,
    const scalar& r,
    const scalar& rb,
    const scalar& Ua,
    const scalar& U,
    const scalar& Ub
)
{

    scalar denominator = (ra - r)*(ra - rb)*(r - rb);

    scalar a1 =
    (
        rb*rb*(Ua - U) + ra*ra*(U - Ub) + r*r*(Ub - Ua)
    ) / denominator;

    scalar a2 =
    (
        rb*(U - Ua) + r*(Ua - Ub) + ra*(Ub - U)
    ) / denominator;

    return a1 + 2.0*a2*r;
}


inline scalar brennerPotential::splineInterpolation
(
    const List<scalar>& x,
    const List<scalar>& y,
    const scalar& xNew
)
{
    scalarField a(4,0.0);
    scalarField b(4,0.0);
    scalarField c(4,0.0);
    scalarField d(4,0.0);

    a[3] = y[0];
    a[2] = 0;

    b[3] = y[1];

    c[3] = y[2];

    d[3] = y[3];

    simpleMatrix<scalar> hmm(11);

//     Info << nl << " y = " << y << ", x = " << x << endl;


    label j = 0;

    // z = 1, y[1]
    hmm[j][0] = 1.0;
    hmm[j][1] = 1.0;
    hmm[j][2] = 0.0;
    hmm[j][3] = 0.0;
    hmm[j][4] = 0.0;
    hmm[j][5] = 0.0;
    hmm[j][6] = 0.0;
    hmm[j][7] = 0.0;
    hmm[j][8] = 0.0;
    hmm[j][9] = 0.0;
    hmm[j][10] = 0.0;

    hmm.source()[j] = y[1] - a[3] - a[2];

    j++;

    // z = 1, y' = y'

    hmm[j][0] = 3.0;
    hmm[j][1] = 2.0;
    hmm[j][2] = 0.0;
    hmm[j][3] = 0.0;
    hmm[j][4] = -1.0;
    hmm[j][5] = 0.0;
    hmm[j][6] = 0.0;
    hmm[j][7] = 0.0;
    hmm[j][8] = 0.0;
    hmm[j][9] = 0.0;
    hmm[j][10] = 0.0;

    hmm.source()[j] = -a[2];

    j++;

    // z = 1, y'' = y''

    hmm[j][0] = 6.0;
    hmm[j][1] = 2.0;
    hmm[j][2] = 0.0;
    hmm[j][3] = -2.0;
    hmm[j][4] = 0.0;
    hmm[j][5] = 0.0;
    hmm[j][6] = 0.0;
    hmm[j][7] = 0.0;
    hmm[j][8] = 0.0;
    hmm[j][9] = 0.0;
    hmm[j][10] = 0.0;
   
    hmm.source()[j] = 0.0;

    j++;

    // z = 1, y[3]

    hmm[j][0] = 0.0;
    hmm[j][1] = 0.0;
    hmm[j][2] = 1.0;
    hmm[j][3] = 1.0;
    hmm[j][4] = 1.0;
    hmm[j][5] = 0.0;
    hmm[j][6] = 0.0;
    hmm[j][7] = 0.0;
    hmm[j][8] = 0.0;
    hmm[j][9] = 0.0;
    hmm[j][10] = 0.0;
   
    hmm.source()[j] = y[2] - b[3];


    j++;

    // z = 1, y' = y'

    hmm[j][0] = 0.0;
    hmm[j][1] = 0.0;
    hmm[j][2] = 3.0;
    hmm[j][3] = 2.0;
    hmm[j][4] = 1.0;
    hmm[j][5] = 0.0;
    hmm[j][6] = 0.0;
    hmm[j][7] = -1.0;
    hmm[j][8] = 0.0;
    hmm[j][9] = 0.0;
    hmm[j][10] = 0.0;
   
    hmm.source()[j] = 0.0;


    j++;

    // z = 1, y'' = y''

    hmm[j][0] = 0.0;
    hmm[j][1] = 0.0;
    hmm[j][2] = 6.0;
    hmm[j][3] = 2.0;
    hmm[j][4] = 0.0;
    hmm[j][5] = 0.0;
    hmm[j][6] = -2.0;
    hmm[j][7] = 0.0;
    hmm[j][8] = 0.0;
    hmm[j][9] = 0.0;
    hmm[j][10] = 0.0;
   
    hmm.source()[j] = 0.0;



    j++;

    // (C) z = 1, 

    hmm[j][0] = 0.0;
    hmm[j][1] = 0.0;
    hmm[j][2] = 0.0;
    hmm[j][3] = 0.0;
    hmm[j][4] = 0.0;
    hmm[j][5] = 1.0;
    hmm[j][6] = 1.0;
    hmm[j][7] = 1.0;
    hmm[j][8] = 0.0;
    hmm[j][9] = 0.0;
    hmm[j][10] = 0.0;
   
    hmm.source()[j] = y[3]-c[3];


    j++;

    // z = 1, y' = y'

    hmm[j][0] = 0.0;
    hmm[j][1] = 0.0;
    hmm[j][2] = 0.0;
    hmm[j][3] = 0.0;
    hmm[j][4] = 0.0;
    hmm[j][5] = 3.0;
    hmm[j][6] = 2.0;
    hmm[j][7] = 1.0;
    hmm[j][8] = 0.0;
    hmm[j][9] = 0.0;
    hmm[j][10] = -1.0;
   
    hmm.source()[j] = 0.0;

    j++;

    // z = 1, y'' = y''

    hmm[j][0] = 0.0;
    hmm[j][1] = 0.0;
    hmm[j][2] = 0.0;
    hmm[j][3] = 0.0;
    hmm[j][4] = 0.0;
    hmm[j][5] = 6.0;
    hmm[j][6] = 2.0;
    hmm[j][7] = 0.0;
    hmm[j][8] = 0.0;
    hmm[j][9] = -2.0;
    hmm[j][10] = 0.0;
   
    hmm.source()[j] = 0.0;



    j++;

    // (D) z = 1, 

    hmm[j][0] = 0.0;
    hmm[j][1] = 0.0;
    hmm[j][2] = 0.0;
    hmm[j][3] = 0.0;
    hmm[j][4] = 0.0;
    hmm[j][5] = 0.0;
    hmm[j][6] = 0.0;
    hmm[j][7] = 0.0;
    hmm[j][8] = 1.0;
    hmm[j][9] = 1.0;
    hmm[j][10] = 1.0;
   
    hmm.source()[j] = y[4]-d[3];

    j++;

    // (D) z = 1, 

    hmm[j][0] = 0.0;
    hmm[j][1] = 0.0;
    hmm[j][2] = 0.0;
    hmm[j][3] = 0.0;
    hmm[j][4] = 0.0;
    hmm[j][5] = 0.0;
    hmm[j][6] = 0.0;
    hmm[j][7] = 0.0;
    hmm[j][8] = 3.0;
    hmm[j][9] = 2.0;
    hmm[j][10] = 1.0;
   
    hmm.source()[j] = 0.0;

    scalarField m(11,0.0);

    m = hmm.LUsolve();

    a[0]=m[0];
    a[1]=m[1];

    b[0]=m[2];
    b[1]=m[3];
    b[2]=m[4];

    c[0]=m[5];
    c[1]=m[6];
    c[2]=m[7];

    d[0]=m[8];
    d[1]=m[9];
    d[2]=m[10];


//     Info << " m = " << m << endl;

    scalar yNew = 0.0;

    if((xNew >= 0) && (xNew <= 1))
    {
        scalar x = 1 - xNew;
        yNew = (a[0]*x*x*x) +(a[1]*x*x) + (a[2]*x) + a[3];
    }

    else if((xNew >= 1) && (xNew <=2))
    {
        scalar x = 2 - xNew;
        yNew = (b[0]*x*x*x) +(b[1]*x*x) + (b[2]*x) + b[3];
    }


    else if((xNew >= 2) && (xNew <= 3))
    {
        scalar x = 3 - xNew;
        yNew = (c[0]*x*x*x) +(c[1]*x*x) + (c[2]*x) + c[3];
    }

    else if((xNew >= 3) && (xNew <= 4))
    {
        scalar x = 4 - xNew;
        yNew = (d[0]*x*x*x) +(d[1]*x*x) + (d[2]*x) + d[3];
    }

    //scalarField &c,

    return yNew;
}


inline scalar brennerPotential::splineInterpolationII
(
    const List<scalar>& x,
    const List<scalar>& y,
    const scalar& xNew
)
{
    scalarField a(4,0.0);
    scalarField b(4,0.0);
    scalarField c(4,0.0);


    a[3] = y[0];
    a[2] = 0;

    b[3] = y[1];

    c[3] = y[2];

    simpleMatrix<scalar> hmm(8);

    label j = 0;

    // z = 1, y[1]
    hmm[j][0] = 1.0;
    hmm[j][1] = 1.0;
    hmm[j][2] = 0.0;
    hmm[j][3] = 0.0;
    hmm[j][4] = 0.0;
    hmm[j][5] = 0.0;
    hmm[j][6] = 0.0;
    hmm[j][7] = 0.0;

    hmm.source()[j] = y[1] - a[3] - a[2];

    j++;

    // z = 1, y' = y'

    hmm[j][0] = 3.0;
    hmm[j][1] = 2.0;
    hmm[j][2] = 0.0;
    hmm[j][3] = 0.0;
    hmm[j][4] = -1.0;
    hmm[j][5] = 0.0;
    hmm[j][6] = 0.0;
    hmm[j][7] = 0.0;

    hmm.source()[j] = -a[2];

    j++;

    // z = 1, y'' = y''

    hmm[j][0] = 6.0;
    hmm[j][1] = 2.0;
    hmm[j][2] = 0.0;
    hmm[j][3] = -2.0;
    hmm[j][4] = 0.0;
    hmm[j][5] = 0.0;
    hmm[j][6] = 0.0;
    hmm[j][7] = 0.0;
   
    hmm.source()[j] = 0.0;

    j++;

    // z = 1, y[3]

    hmm[j][0] = 0.0;
    hmm[j][1] = 0.0;
    hmm[j][2] = 1.0;
    hmm[j][3] = 1.0;
    hmm[j][4] = 1.0;
    hmm[j][5] = 0.0;
    hmm[j][6] = 0.0;
    hmm[j][7] = 0.0;
   
    hmm.source()[j] = y[2] - b[3];


    j++;

    // z = 1, y' = y'

    hmm[j][0] = 0.0;
    hmm[j][1] = 0.0;
    hmm[j][2] = 3.0;
    hmm[j][3] = 2.0;
    hmm[j][4] = 1.0;
    hmm[j][5] = 0.0;
    hmm[j][6] = 0.0;
    hmm[j][7] = -1.0;
   
    hmm.source()[j] = 0.0;


    j++;

    // z = 1, y'' = y''

    hmm[j][0] = 0.0;
    hmm[j][1] = 0.0;
    hmm[j][2] = 6.0;
    hmm[j][3] = 2.0;
    hmm[j][4] = 0.0;
    hmm[j][5] = 0.0;
    hmm[j][6] = -2.0;
    hmm[j][7] = 0.0;

   
    hmm.source()[j] = 0.0;

    j++;


    hmm[j][0] = 0.0;
    hmm[j][1] = 0.0;
    hmm[j][2] = 0.0;
    hmm[j][3] = 0.0;
    hmm[j][4] = 0.0;
    hmm[j][5] = 1.0;
    hmm[j][6] = 1.0;
    hmm[j][7] = 1.0;

   
    hmm.source()[j] = y[3]-c[3];

    j++;

    hmm[j][0] = 0.0;
    hmm[j][1] = 0.0;
    hmm[j][2] = 0.0;
    hmm[j][3] = 0.0;
    hmm[j][4] = 0.0;
    hmm[j][5] = 3.0;
    hmm[j][6] = 2.0;
    hmm[j][7] = 1.0;
   
    hmm.source()[j] = 0.0;


    scalarField m(8,0.0);

    m = hmm.LUsolve();

    a[0]=m[0];
    a[1]=m[1];

    b[0]=m[2];
    b[1]=m[3];
    b[2]=m[4];

    c[0]=m[5];
    c[1]=m[6];
    c[2]=m[7];



    scalar yNew = 0.0;

    if((xNew >= 0) && (xNew <= 1))
    {
        scalar x = 1 - xNew;
        yNew = (a[0]*x*x*x) +(a[1]*x*x) + (a[2]*x) + a[3];
    }

    else if((xNew >= 1) && (xNew <= 2))
    {
        scalar x = 2 - xNew;
        yNew = (b[0]*x*x*x) +(b[1]*x*x) + (b[2]*x) + b[3];
    }

    else if((xNew >= 2) && (xNew <= 3))
    {
        scalar x = 3 - xNew;
        yNew = (c[0]*x*x*x) +(c[1]*x*x) + (c[2]*x) + c[3];
    }

    //scalarField &c,

    return yNew;
}

// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
