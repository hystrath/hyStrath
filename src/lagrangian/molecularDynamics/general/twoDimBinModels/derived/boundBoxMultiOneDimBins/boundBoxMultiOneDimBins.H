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

    Create a set of bins within a cuboid domain from start point to end point.
    Each bin is further divided in a perpendicular direction d1 into another 
    set of bins,where a one dimensional profile is measured. The result are N
    number of one dimensional profiles.
    

SourceFiles
    boundBoxMultiOneDimBinsI.H
    boundBoxMultiOneDimBins.C
    boundBoxMultiOneDimBinsIO.C

\*---------------------------------------------------------------------------*/

#ifndef boundBoxMultiOneDimBins_H
#define boundBoxMultiOneDimBins_H


#include "twoDimBinModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class boundBoxMultiOneDimBins Declaration
\*---------------------------------------------------------------------------*/

class boundBoxMultiOneDimBins
:
    public twoDimBinModel
{

private:

    // Private data

        dictionary propsDict_;

        
        
        // outer bound box
        
        dictionary outerDict_;
        
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
       
        
//         scalar rSEMag_;
        
        boundBox outerBB_;
        
        List<vector> rBins_;        

        
        // inner bound boxes
        
        dictionary innerDict_;
        
       
        label nBinsY_;
        vector unitVectorY_;

        vector startPointY_;
        vector endPointY_;        
        
        scalar binWidthY_;
        
        List<boundBox> innerBBs_;
        

        
        
    // Private Member Functions

        void checkBoundBox
        (
            boundBox& b,
            const vector& startPoint,
            const vector& endPoint
        );

public:

    //- Runtime type information
    TypeName("boundBoxMultiOneDimBins");

    // Static data members

    // Constructors

        //- construct from components
        boundBoxMultiOneDimBins
        (
            const polyMesh& mesh,
            const dictionary& dict
        );

    // Selectors


    // Destructor

        ~boundBoxMultiOneDimBins();

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
