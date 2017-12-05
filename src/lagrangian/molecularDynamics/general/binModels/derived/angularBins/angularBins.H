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
    angularBinsI.H
    angularBins.C
    angularBinsIO.C

\*---------------------------------------------------------------------------*/

#ifndef angularBins_H
#define angularBins_H


#include "binModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class angularBins Declaration
\*---------------------------------------------------------------------------*/

class angularBins
:
    public binModel
{

private:

    // Private data

        dictionary propsDict_;

//         vector rotationalAxis_;

        // midpoint between rotational start and end points
        vector m_;
        
        // define two reference vectors orthogonal to each other and the rotational axis
        vector rVx_;
        vector rVy_;
        
        // radius 
        scalar Rin_;
        scalar Rout_;

        // angle with respect to the reference vectors (important to define clockwise)
        scalar thetaStart_;
        scalar thetaEnd_;      

//         vector startPoint_;
//         vector endPoint_;
//         vector unitVector_;
//         scalar rSEMag_;


        label nBins_;

        // in radians
        scalar binWidth_;
        scalar area_;

        bool counterClockWise_;


public:

    //- Runtime type information
    TypeName("angularBins");

    // Static data members

    // Constructors

        //- construct from components
        angularBins
        (
            const polyMesh& mesh,
            const dictionary& dict
        );

    // Selectors


    // Destructor

        ~angularBins();

    // Member Functions

        // Access

        // Check

        // Edit

            label isPointWithinBin
            (
                const vector& rI,
                const label& cellI
            );

            scalarField binPositions();

            vectorField bins();

            const label& nBins() const;

            scalar binVolume(const label& n);

            

        // Write





    // Member Operators

//         void operator=(const angularBins&);

    // Friend Functions

    // Friend Operators

    // IOstream Operators

//         friend Ostream& operator<<(Ostream&, const angularBins&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
