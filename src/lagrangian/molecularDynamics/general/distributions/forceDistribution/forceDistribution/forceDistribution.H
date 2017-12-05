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
    A linear force distribution. Can be read from file or constructed
    appropriately.


Description

SourceFiles
    forceDistributionI.H
    forceDistribution.C
    forceDistributionIO.C

\*---------------------------------------------------------------------------*/

#ifndef forceDistribution_H
#define forceDistribution_H

#include "moleculeCloud.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class forceDistribution Declaration
\*---------------------------------------------------------------------------*/

class forceDistribution
{
    // Private data

        word name_;

        word fileName_;

        vector startPoint_;

        vector endPoint_;

        vector unitVector_;

        label noOfBins_;

        vectorField forces_;
        scalarField energies_;
        scalarField virials_;

        scalarField mols_;

        vectorField radii_;
        scalarField magRadii_;

        scalarField magForces_;

        scalar binWidth_;

    // Private Member Functions

        void setRadii();

public:

    // Constructors

        //- Construct null
        forceDistribution();

        //- Construct from name
        forceDistribution
        (
            const word& name
//             const word& fileName
        );

        //- Construct from components
        forceDistribution
        (
            const word& name,
//             const word& fileName,
            const vector& startPoint,
            const vector& endPoint,
            const label& noOfBins
        );

        forceDistribution
        (
            const dictionary& dict
        );

    // Destructor

        ~forceDistribution();

    // Member Functions

        void setForceDistr
        (
            const word& name,
//             const word& fileName,
            const vector& startPoint,
            const vector& endPoint,
            const label& noOfBins
        );

        void setForceDistr
        (
            const dictionary& dict
        );

        void setMagForceDistr();

        void scaleForceDistribution(const scalar& value);

        void scaleForceDistribution(const scalarField& values);

        void scaleSampledDistribution();

        void setForceDistribution
        (
            const List< Pair<scalar> >&
        );

        void setPEDistribution
        (
            const List< Pair<scalar> >&
        );

        void setVirialDistribution
        (
            const List< Pair<scalar> >&
        );

        scalar readBinWidth();

        void addToDistribution
        (
            const vector& r, 
            const vector& force,
            const scalar& energy
        );

        void addToDistribution
        (
            const vector& r, 
            const vector& force,
            const scalar& energy,
            const scalar& virial
        );

        bool isWithinDistributionRange(const vector& r);

        // Access

            List< Pair<scalar> > fMagDistr();
            List< Pair<scalar> > pEpairDistr();
            List< Pair<scalar> > virialPairDistr();

            inline const word& name() const;
            inline const word& filename() const;
            inline const label& noOfBins() const;
            inline const scalar& binWidth() const;

            inline const vector& startPoint() const;
            inline const vector& endPoint() const;
            inline const vector& unitVector() const;

            inline const vectorField& radii() const;
            inline const scalarField& magRadii() const;

            inline const vectorField& forces() const;
            inline vectorField& forces();

            inline const scalarField& magForces() const;

            inline const scalarField& energies() const;
            inline scalarField& energies();


            // returns bin values for force, PE and virial
            scalar binValue(const scalar& r) const;
            scalar binValuePE(const scalar& rD) const;
            scalar binValueVirial(const scalar& rD) const;

            // returns linear interpolated value for force, PE and virial
            scalar linearInterp(const scalar& r) const;
            scalar linearInterpPE(const scalar& r) const;
            scalar linearInterpVirial(const scalar& r) const;



        // Check

        // Edit

        // Write

            void clear();

            void read(const Time& runTime);

            void write(const Time& runTime);

            void writeLinearForceDistribution(const Time& runTime);

            void writeTimeData
            (
                const fileName& pathName,
                const word& nameFile,
                const scalarField& xData,
                const scalarField& yData
            );

            //- vector field
            void writeTimeData
            (
                const fileName& pathName,
                const word& nameFile,
                const scalarField& xData,
                const vectorField& yData
            );


//             void setInitialPBCDistribtion
//             (
//                 const moleculeCloud& mD
//             );

//             void setForcesFromPBC
//             (
//                 const moleculeCloud& mD,
//                 const label& molId,
//                 const scalar& rMin
//             );


    // Member Operators

//         void operator=(const forceDistribution&);

    // Friend Functions

    // Friend Operators

    // IOstream Operators

//         friend Ostream& operator<<(Ostream&, const forceDistribution&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "forceDistributionI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
