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
    simpleBinsI.H
    simpleBins.C
    simpleBinsIO.C

\*---------------------------------------------------------------------------*/

#ifndef simpleBins_H
#define simpleBins_H


#include "binModel.H"
#include "boundedBox.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class simpleBins Declaration
\*---------------------------------------------------------------------------*/

class simpleBins
:
    public binModel
{

private:

    // Private data

        dictionary propsDict_;

        vector unitVector_;

        boundedBox bb_;

        vector startPoint_;
        vector endPoint_;

        scalar rSEMag_;


        label nBins_;
        scalar binWidth_;
        scalar area_;




public:

    //- Runtime type information
    TypeName("simpleBins");

    // Static data members

    // Constructors

        //- construct from components
        simpleBins
        (
            const polyMesh& mesh,
            const dictionary& dict
        );

    // Selectors


    // Destructor

        ~simpleBins();

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

//         void operator=(const simpleBins&);

    // Friend Functions

    // Friend Operators

    // IOstream Operators

//         friend Ostream& operator<<(Ostream&, const simpleBins&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
