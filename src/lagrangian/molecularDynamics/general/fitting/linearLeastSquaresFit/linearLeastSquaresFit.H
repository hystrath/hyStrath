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


Description
    Taken from Numerical Recipes in C - The art of scientific computing 2nd edition
    William H. Press
    Pg 665

SourceFiles
    linearLeastSquaresFitI.H
    linearLeastSquaresFit.C
    linearLeastSquaresFitIO.C

\*---------------------------------------------------------------------------*/

#ifndef linearLeastSquaresFit_H
#define linearLeastSquaresFit_H

// #include "distribution.H"
#include "scalarField.H"
#include "polyMesh.H"
#include "IFstream.H"
#include "OFstream.H"
#include "Time.H"
// #include <math.h>
// #include "nrutil.h"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class linearLeastSquaresFit Declaration
\*---------------------------------------------------------------------------*/

class linearLeastSquaresFit
{
    // Private data


    // Private Member Functions

//         void setRadius();

public:

    // Static data members

    // Constructors

        //- Construct null
//         linearLeastSquaresFit();

        //- Construct from xs and ys
        linearLeastSquaresFit
        (
            const scalarField& x,
            const scalarField& y,
//             scalarField sig,
//             label mwt,
            scalar& a,
            scalar& b
//             scalar& siga,
//             scalar& sigb,
//             scalar& chi2,
//             scalar& q
        );


    // Selectors

    // Destructor

        ~linearLeastSquaresFit();

    // Member Functions



    // Member Operators

//         void operator=(const linearLeastSquaresFit&);

    // Friend Functions

    // Friend Operators

    // IOstream Operators

//         friend Ostream& operator<<(Ostream&, const linearLeastSquaresFit&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
