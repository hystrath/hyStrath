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

    no-time counter with reactions

SourceFiles
    simplifiedBernoulliTrialsI.H
    simplifiedBernoulliTrials.C
    simplifiedBernoulliTrialsIO.C

\*---------------------------------------------------------------------------*/

#ifndef simplifiedBernoulliTrialsTAS_H
#define simplifiedBernoulliTrialsTAS_H


#include "collisionPartnerSelection.H"
#include "vectorField.H"
#include "dsmcCloud.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class simplifiedBernoulliTrialsTAS Declaration
\*---------------------------------------------------------------------------*/

class simplifiedBernoulliTrialsTAS
:
    public collisionPartnerSelection
{

private:

    // Private data

        dictionary propsDict_;

        // number of parcels per sub-cell
        label n_;

        // adaptive time-steps after each measurement
        label nSteps_;
//         scalar subCellVolume_;

        label counter_;

        // measured parcels per volume
        scalarField rhoN_;

//         scalarField subCellVolumes_;

        vectorField startPoints_;

        // nx, ny nz slices in x y z directions 
        labelListList nSlices_;

        // no of actual sub-cells 
        labelList nSubCells_;



        // x, y, z components of vector is deltax, deltay and deltaz
        vectorField binWidths_;

		    // Private Member Functions

            pointField cellPoints(const label& cell);


            void measureLocalDensity();

            // check for skewness etc
            void checkSubCelling();


public:

    //- Runtime type information
    TypeName("simplifiedBernoulliTrialsTAS");

    // Static data members

    // Constructors

        //- construct from components
        simplifiedBernoulliTrialsTAS
        (
            const polyMesh& mesh,
            dsmcCloud& cloud,
            const dictionary& dict
        );

    // Selectors


    // Destructor

        ~simplifiedBernoulliTrialsTAS();

    // Member Functions

        // Access

        // Check

        // Edit

            void initialConfiguration();

            void collide();

        // Write





    // Member Operators

//         void operator=(const simplifiedBernoulliTrials&);

    // Friend Functions

    // Friend Operators

    // IOstream Operators

//         friend Ostream& operator<<(Ostream&, const simplifiedBernoulliTrials&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
