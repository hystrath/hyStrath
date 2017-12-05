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

SourceFiles
    leastSquaresFitI.H
    leastSquaresFit.C
    leastSquaresFitIO.C

\*---------------------------------------------------------------------------*/

#ifndef leastSquaresFit_H
#define leastSquaresFit_H

// #include "distribution.H"
#include "scalarField.H"
#include "polyMesh.H"
#include "IFstream.H"
#include "OFstream.H"
#include "Time.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class leastSquaresFit Declaration
\*---------------------------------------------------------------------------*/

class leastSquaresFit
{
    // Private data

        word name_;

        bool write_;

        scalarField xs_;

        scalarField ys_;

        scalar gradient_;

        scalar yIntercept_;


    // Private Member Functions

//         void setRadius();

public:

    // Static data members

    // Constructors

        //- Construct null
        leastSquaresFit();

        //- Construct from xs and ys
        leastSquaresFit
        (
            const word& name,
            const bool& write,
            const scalarField& xs,
            const scalarField& ys
        );


    // Selectors

    // Destructor

        ~leastSquaresFit();

    // Member Functions

        void setInitialData
        (
            const word& name,
            const bool& write,
            const label nBins
        );

        void setFitParameters();

        // Access

            inline const scalar& gradient() const;

            inline const scalar& yIntercept() const;

            inline const scalarField& xs() const;

            inline const scalarField& ys() const;

            inline scalarField& xs();

            inline scalarField& ys();
        // Check

        // Edit

        // Write

            void writeField(const Time& runTime);

            void writeTimeData
            (
                const fileName& pathName,
                const word& nameFile,
                const scalarField& xData,
                const scalarField& yData
            );


    // Member Operators

//         void operator=(const leastSquaresFit&);

    // Friend Functions

    // Friend Operators

    // IOstream Operators

//         friend Ostream& operator<<(Ostream&, const leastSquaresFit&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "leastSquaresFitI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
