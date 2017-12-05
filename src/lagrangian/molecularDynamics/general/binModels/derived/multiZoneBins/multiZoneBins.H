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
    multiZoneBinsI.H
    multiZoneBins.C
    multiZoneBinsIO.C

\*---------------------------------------------------------------------------*/

#ifndef multiZoneBins_H
#define multiZoneBins_H


#include "binModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class multiZoneBins Declaration
\*---------------------------------------------------------------------------*/

class multiZoneBins
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

        List<label> regionIds_;
        List<scalar> volumes_;

        List<scalar> binWidths_;

        labelField cellRegionAddressing_;


public:

    //- Runtime type information
    TypeName("multiZoneBins");

    // Static data members

    // Constructors

        //- construct from components
        multiZoneBins
        (
            const polyMesh& mesh,
            const dictionary& dict
        );

    // Selectors


    // Destructor

        ~multiZoneBins();

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

//         void operator=(const multiZoneBins&);

    // Friend Functions

    // Friend Operators

    // IOstream Operators

//         friend Ostream& operator<<(Ostream&, const multiZoneBins&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
