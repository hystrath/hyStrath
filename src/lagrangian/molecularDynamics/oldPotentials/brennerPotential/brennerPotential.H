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
    OBSOLETE - PROPER IMPLEMENTATION NOT FOUND


SourceFiles
    brennerPotentialI.H
    brennerPotential.C
    brennerPotentialIO.C

\*---------------------------------------------------------------------------*/

#ifndef brennerPotential_H
#define brennerPotential_H

#include "polyMesh.H"
#include "potential.H"
#include "dictionary.H"
// #include "writeTimeData.H"
#include "mathematicalConstants.H"
#include "reducedUnits.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class brennerPotential Declaration
\*---------------------------------------------------------------------------*/

class brennerPotential
{
    // Private data

        const polyMesh& mesh_;
        const reducedUnits& rU_;
//         const dictionary& potentialProperties_;

        // potential constants 
        scalar sIJ_;
        scalar cIJ_;
        scalar betaIJ_;
        scalar rIJ0_;
        scalar rIJ1_;
        scalar rIJ2_;
        scalar deltaI_;
        scalar alphaIJK_;
        scalar a0_;
        scalar c0_;
        scalar d0_;

        scalar dr_;
        scalar rMin_;
        scalar rCut_;
        scalar constantRepulsive_;
        scalar constantAttractive_;

        scalar constantBrenner_;

    // Private Member Functions
        

public:

    // Static data members

    // Constructors

        //- Construct null
        brennerPotential
        (
            const polyMesh& mesh,
            const reducedUnits& redUnits
        );



    // Selectors

    // Destructor

        ~brennerPotential();

    // Member Functions

        void setPotential();

        // Access

//             inline List<label> molIds();
//             inline labelList& molIds();

        inline const scalar& rCutMax() const;

        inline scalar Urep(const scalar& r);
        inline scalar Uatt(const scalar& r);

        inline scalar fRep(const scalar& r);
        inline scalar fAtt(const scalar& r);


        inline scalar forceRepulsive(const scalar& r);
        inline scalar forceAttractive(const scalar& r);


//         inline scalar force(const scalar& r, const scalar& Bij);
//         inline scalar potentialEnergy(const scalar& r, const scalar& Bij);

        inline scalar fIJ(const scalar& r);

        inline scalar F(const scalar& x);

        inline scalar G(const scalar& cosThetaIJK);

        inline scalar numericalDifferentiation
        (
            const scalar& ra,
            const scalar& r,
            const scalar& rb,
            const scalar& Ua,
            const scalar& U,
            const scalar& Ub
        );

        inline scalar splineInterpolation
        (
            const List<scalar>& x,
            const List<scalar>& y,
            const scalar& xNew
        );

        inline scalar splineInterpolationII
        (
            const List<scalar>& x,
            const List<scalar>& y,
            const scalar& xNew
        );

        // Check

        // Edit

        // Write
            void outputProperties();



    // Member Operators

//         void operator=(const brennerPotential&);

    // Friend Functions

    // Friend Operators

    // IOstream Operators

//         friend Ostream& operator<<(Ostream&, const brennerPotential&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "brennerPotentialI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
