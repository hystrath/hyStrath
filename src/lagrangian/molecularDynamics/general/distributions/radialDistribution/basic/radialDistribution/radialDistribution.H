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
    radialDistributionI.H
    radialDistribution.C
    radialDistributionIO.C

\*---------------------------------------------------------------------------*/

#ifndef radialDistribution_H
#define radialDistribution_H

// #include "distribution.H"
// #include "moleculeCloud.H"
#include "IFstream.H"
#include "OFstream.H"
#include "polyMesh.H"
#include "vector.H"
#include "volFields.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class radialDistribution Declaration
\*---------------------------------------------------------------------------*/

class radialDistribution
// :
//     public distribution
{
    // Private data

        word name_;

        word fileName_;

        scalar rMax_;

        label noOfBins_;

        scalar binWidth_;

        scalarField g_;

        scalarField mols_;

        scalarField radius_;





    // Private Member Functions

        void setRadius();

public:

    // Static data members

    // Constructors

        //- Construct null
        radialDistribution();

        //- Construct from name
        //- (used for reading an rdf from input file)
        radialDistribution(const word& name);

        //- Construct from name, rMax noOfBins
        //- (used for sampling radial distribution)
        radialDistribution
        (
            const word& name,
            const scalar& rMax,
            const label& noOfBins
        );

        //- construct from components
        //- initial conditions or for cases in which rdf remains fixed
        radialDistribution
        (
            const word& name,
            const List< scalar >& g,
            const List< scalar >& radius
        );

        //- construct from components: pair lists
        radialDistribution
        (
            const word& name,
            const List< Pair<scalar> >& rdf
        );


    // Selectors

    // Destructor

        ~radialDistribution();

    // Member Functions

        void setRDF
        (
            const scalar& rMax,
            const label& noOfBins
        );

        void setRdfName
        (
            const word& name
        );

        void setRdf
        (
            const word& name,
            const List< scalar >& g,
            const List< scalar >& radius
        );

        void setRDF(const List< Pair<scalar> >&);

        void setRdf
        (
            const scalar& rhoAv,
            const scalar& nMols
        );

        // returns binWidth from a read distribution
        scalar readBinWidth();

        void addToDistribution(const scalar& r);

        void clearRadialDistribution();

        // Access

            // returns linear interpolated g
            scalar gLinear(const scalar& r) const;

            // returns bin value for g
            scalar gBin(const scalar& r) const;

            List< Pair<scalar> > radialDistrib();

            inline const scalarField& g() const;
            inline scalarField& g();

            inline const scalarField& r() const;

            inline const word& name() const;

            inline const scalar& rMax() const;

            inline const scalar& binWidth() const;

            inline const scalarField& mols() const;

            void clearCollector();
        // Check

        // Edit

        // Write

            void writeRDF(const Time& runTime);

            void readRDF(const Time& runTime);

            void writeTimeData
            (
                const fileName& pathName,
                const word& nameFile,
                const scalarField& xData,
                const scalarField& yData
            );


    // Member Operators

//         void operator=(const radialDistribution&);

    // Friend Functions

    // Friend Operators

    // IOstream Operators

//         friend Ostream& operator<<(Ostream&, const radialDistribution&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "radialDistributionI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
