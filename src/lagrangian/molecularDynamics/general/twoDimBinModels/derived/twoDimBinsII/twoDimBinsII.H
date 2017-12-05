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
    startPoint and endPoint define the boundbox.
    unit vectors need to be defined from startPoint:
    
    
SourceFiles
    twoDimBinsIII.H
    twoDimBinsII.C
    twoDimBinsIIIO.C

\*---------------------------------------------------------------------------*/

#ifndef twoDimBinsII_H
#define twoDimBinsII_H


#include "twoDimBinModel.H"
#include "boundedBox.H" 

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class twoDimBinsII Declaration
\*---------------------------------------------------------------------------*/

class twoDimBinsII
:
    public twoDimBinModel
{

private:

    // Private data

        dictionary propsDict_;

        vector startPoint_;
        vector endPoint_;
        
        boundedBox box_;
                
        vector unitVectorX_;
        vector unitVectorY_;
        vector unitVectorZ_; // thickness
        
        label nBinsX_;
        label nBinsY_;

        scalar lengthX_;
        scalar lengthY_;
        scalar lengthZ_; // thickness
        
        scalar binWidthX_;
        scalar binWidthY_;
        
        
    // Private Member Functions

        void checkBoundBox
        (
            boundBox& b,
            const vector& startPoint,
            const vector& endPoint
        );
        
        
public:

    //- Runtime type information
    TypeName("twoDimBinsII");

    // Static data members

    // Constructors

        //- construct from components
        twoDimBinsII
        (
            const polyMesh& mesh,
            const dictionary& dict
        );

    // Selectors


    // Destructor

        ~twoDimBinsII();

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
