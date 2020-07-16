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

SourceFiles
    collisionDistributionI.H
    collisionDistribution.C
    collisionDistributionIO.C

\*---------------------------------------------------------------------------*/

#ifndef collisionDistribution_H
#define collisionDistribution_H

#include "distribution.H"
// #include "moleculeCloud.H"
#include "IFstream.H"
#include "OFstream.H"
#include "polyMesh.H"
#include "vector.H"
#include "volFields.H"
#include "writeTimeData.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class collisionDistribution Declaration
\*---------------------------------------------------------------------------*/

class collisionDistribution
:
    public distribution
{
    // Private data

        word name_;

        scalar rMax_;

        scalarField p_;

        scalarField radius_;



    // Private Member Functions

        void setRadius();

public:

    // Static data members

    // Constructors

        //- Construct null
        collisionDistribution();

        //- Construct from name
        //- (used for reading an rdf from input file)
        collisionDistribution(const word&);

        //- Construct from binWidth, rMax and name
        collisionDistribution
        (
            const scalar&,
            const scalar&,
            const word&
        );


    // Selectors

    // Destructor

        ~collisionDistribution();

    // Member Functions

        void setDistribution();

        void clearRadialDistribution();

        // Access

            // returns linear interpolated g
            scalar gLinear(const scalar& r) const;

            // returns bin value for g
            scalar gBin(const scalar& r) const;

            List< Pair<scalar> > distrib();

            inline const List<scalar>& p() const;
            inline scalarField& p();

            inline const scalarField& r() const;

            inline const word& name() const;

            inline const scalar& rMax() const;

            void addPoint(const scalar& r);

        // Check

        // Edit

        // Write

            void writeDistribution(const Time& runTime);
//
//             void writeTimeData
//             (
//                 const fileName& pathName,
//                 const word& nameFile,
//                 const scalarField& xData,
//                 const scalarField& yData
//             );
    // Member Operators

//         void operator=(const collisionDistribution&);

    // Friend Functions

    // Friend Operators

    // IOstream Operators

//         friend Ostream& operator<<(Ostream&, const collisionDistribution&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "collisionDistributionI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
