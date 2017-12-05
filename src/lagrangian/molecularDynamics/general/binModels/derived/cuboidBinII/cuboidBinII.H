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
    cuboidBinIII.H
    cuboidBinII.C
    cuboidBinIIIO.C

\*---------------------------------------------------------------------------*/

#ifndef cuboidBinII_H
#define cuboidBinII_H

#include "binModel.H"
#include "boundBox.H" 

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class cuboidBinII Declaration
\*---------------------------------------------------------------------------*/

class cuboidBinII
:
    public binModel
{

private:

    // Private data

        dictionary propsDict_;

        vector startPoint_;
        vector endPoint_;
        vector unitVector_;

        label nBins_;
        scalar binWidth_;

        scalar d1_;
        scalar d2_;
        vector n1_;
        vector n2_;

        scalar area_;

        scalar rSEMag_;
        vector midpoint_;
        
        boundBox box_;
        
        
        // Private Member Functions
            
        void checkBoundBox
        (
            boundBox& b,
            const vector& startPoint,
            const vector& endPoint
        );     
        
public:

    //- Runtime type information
    TypeName("cuboidBinII");

    // Static data members

    // Constructors

        //- construct from components
        cuboidBinII
        (
            const polyMesh& mesh,
            const dictionary& dict
        );

    // Selectors


    // Destructor

        ~cuboidBinII();

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

//         void operator=(const cuboidBinII&);

    // Friend Functions

    // Friend Operators

    // IOstream Operators

//         friend Ostream& operator<<(Ostream&, const cuboidBinII&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
