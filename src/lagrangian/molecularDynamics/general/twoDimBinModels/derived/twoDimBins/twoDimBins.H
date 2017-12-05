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
    twoDimBinsI.H
    twoDimBins.C
    twoDimBinsIO.C

\*---------------------------------------------------------------------------*/

#ifndef twoDimBins_H
#define twoDimBins_H


#include "twoDimBinModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class twoDimBins Declaration
\*---------------------------------------------------------------------------*/

class twoDimBins
:
    public twoDimBinModel
{

private:

    // Private data

        dictionary propsDict_;

        vector startPoint_;
        vector endPoint_;
        scalar h_;
        vector unitVector_;

        vector unitVectorX_;
        vector unitVectorY_;
        
        label nBinsX_;
        label nBinsY_;

        scalar dX_;
        scalar dY_;
        
        scalar binWidthX_;
        scalar binWidthY_;
        
        vector startPointX_;
        vector startPointY_;
        
        vector endPointX_;
        vector endPointY_;
        
        
        
    // Private Member Functions


public:

    //- Runtime type information
    TypeName("twoDimBins");

    // Static data members

    // Constructors

        //- construct from components
        twoDimBins
        (
            const polyMesh& mesh,
            const dictionary& dict
        );

    // Selectors


    // Destructor

        ~twoDimBins();

    // Member Functions

        // Access

        // Check

        // Edit

            // 0 length, 1 radius, 2 angle
            List<label> isPointWithinBin
            (
                const vector& rI,
                const label& cellI
            );

            scalarField binPositionsX();
            scalarField binPositionsY();
//             scalarField binPositionsZ();

            vector position(/*const vector&, */const scalar&, const scalar&);

//             Field <Field<vectorField > > bins();
//             vectorField bins();

            // 0 length, 1 radius, 2 angle
            List<label> nBins();

            scalar binVolume(const label& n);

            void write
            (
                const fileName& path,
                const word& name
            );

            const vector& refVector() const;

        // Write





    // Member Operators

    // Friend Functions

    // Friend Operators

    // IOstream Operators


};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
