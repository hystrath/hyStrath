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
    twoDimForceDistributionI.H
    twoDimForceDistribution.C
    twoDimForceDistributionIO.C

\*---------------------------------------------------------------------------*/

#ifndef twoDimForceDistribution_H
#define twoDimForceDistribution_H


// #include "distribution.H"
#include "moleculeCloud.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class twoDimForceDistribution Declaration
\*---------------------------------------------------------------------------*/

class twoDimForceDistribution
{
    // Private data

        word name_;
        word fileName_;

        vector startPoint_;

        scalar lengthX_;
        scalar lengthY_;

        vector unitVectorX_;
        vector unitVectorY_;

        label nBinsX_;
        label nBinsY_;

        scalar binWidthX_;
        scalar binWidthY_;

//         scalarField magRadiiX_;
//         scalarField magRadiiY_;

        List< scalarField > mols_;
        List< vectorField > forces_;
        List< scalarField > energies_;

    // Private Member Functions

//         void setRadii();

public:

    // Constructors

        //- Construct null
        twoDimForceDistribution();

        //- Construct from dictionary
        twoDimForceDistribution
        (
            const dictionary& dict
        );

    // Destructor

        ~twoDimForceDistribution();

    // Member Functions

        void scaleForceDistribution(const scalar& value);

        void scaleForceDistribution(const List<scalarField>& values);

        void scaleSampledDistribution();

        void addToDistribution
        (
            const vector& r, 
            const vector& force,
            const scalar& energy
        );

        bool isWithinDistributionRange(const vector& r);

        // Access

//             List< Pair<scalar> > fMag();

            inline const word& name() const;
            inline const word& filename() const;
            inline const vector& startPoint() const;

            inline const vector& unitVectorX() const;
            inline const vector& unitVectorY() const;

            inline const scalar& binWidthX() const;
            inline const scalar& binWidthY() const;

            inline const label& nBinsX() const;
            inline const label& nBinsY() const;

//             inline const vectorField& forces() const;
//             inline vectorField& forces();
// 
//             inline const scalarField& energies() const;
//             inline scalarField& energies();

        // Check

        // Edit

        // Write

            void clear();

            void write
            (
                const Time& runTime,
                const polyMesh& mesh
            );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "twoDimForceDistributionI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
