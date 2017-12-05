/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2007 OpenCFD Ltd.
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
    integral

Description

	Stores all the information for the integral models (state and fluxes),
    defined within the integralDict, and selects & builds the models automatically.

\*---------------------------------------------------------------------------*/

#ifndef integral_H
#define integral_H

#include "vector.H"
#include "List.H"
// #include "IOdictionary.H"
// #include "dictionary.H"
#include "autoPtr.H"
#include "polyMesh.H"
#include "volFields.H"
// #include "Time.H"
// #include "OFstream.H"
// #include "IFstream.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class integral Declaration
\*---------------------------------------------------------------------------*/

class integral
{

    // Private data

        scalar area_;


public:

    // Constructors

        integral
        (
            const scalar& binWidth,
            const scalarField& xData,
            const scalarField& yData,
            const word& integrationType
        );



        ~integral();


    // Member Functions

        const scalar& area() const;
    // Ostream Operator

//         friend Ostream& operator<<(Ostream&, const integral&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
