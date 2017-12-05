/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2007 OpenCFD Ltd.
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
    polyTwoDimDropletDensity

Description
    two dimensional measurement of the contour lines of a droplet

    Bins divided into height and radial components
    - one bin is a circular bin with a rectangular cross-section

SourceFiles
    polyTwoDimDropletDensityI.H
    polyTwoDimDropletDensity.C
    polyTwoDimDropletDensityIO.C

\*---------------------------------------------------------------------------*/

#ifndef polyTwoDimDropletDensity_H
#define polyTwoDimDropletDensity_H

#include "polyField.H"
#include "vector.H"
#include "polyMesh.H"
#include "volFields.H"
#include "polyMoleculeCloud.H"
#include "distribution.H"
#include "boundedBox.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class polyTwoDimDropletDensity Declaration
\*---------------------------------------------------------------------------*/

class polyTwoDimDropletDensity
:
	public polyField
{

private:

    // Private data

        //- subDictionary containing the properties
        dictionary propsDict_;

        propertyField fields_;

        word fieldName_;

       
        // total number of species in droplet
        List<label> dropletMolIds_;
        
        // species you want to measure
        List<label> molIds_;
        
        
        // bottom centre point of droplet
        vector startPoint_;
        
        // top centre point of droplet
        vector endPoint_;

        scalar h_;

        vector unitVectorX_;

        scalar radius_;

        // no of bins in the height dir h
        label nBinsX_;
        
        // no of bins in the radial dir R (may vary depending on algorithm)
        label nBinsY_;

        scalar binWidthX_;
        scalar binWidthY_;
        
        vector oldCenterOfMass_;
        
        
        Field<scalarField> mols_;
        Field<scalarField> mass_;

        Field<scalarField> rhoN_;
        Field<scalarField> rhoM_;
        
        scalarField magRadii_;
        scalarField binWidths_;
        scalarField volume_;

        scalar avVolume_;

        // for searching
        scalar minBinWidth_;
        labelField n_;

        scalar nAvTimeSteps_;
        
        boundedBox box_;
        scalar x_;
        scalar y_;
        scalar z_;
        
        bool resetAtOutput_;
        
    // Private Member Functions

        label findBin(const scalar& r);


public:

    //- Runtime type information
    TypeName("polyTwoDimDropletDensity");

    // Constructors

        //- Construct from components
        polyTwoDimDropletDensity
        (
            Time& t,
            const polyMesh& mesh,
            polyMoleculeCloud& molCloud,
            const dictionary& dict
        );

    // Destructor

        ~polyTwoDimDropletDensity();


    // Member Functions

        //- create initial field
        void createField();

        //- calculate field
        void calculateField();

        void measureDuringForceComputation
        (
            polyMolecule* molI,
            polyMolecule* molJ
        );

        void measureDuringForceComputationSite
        (
            polyMolecule* molI,
            polyMolecule* molJ,
            label sI,
            label sJ
        );

        //- return field
//         const volVectorField& velocityField() const;
        const propertyField& fields() const;

        //- write field
        void writeField();

//         void updateProperties(const dictionary&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
