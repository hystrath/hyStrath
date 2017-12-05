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
    polyDensityMOB2DRadialCOM

Description


SourceFiles
    polyDensityMOB2DRadialCOMI.H
    polyDensityMOB2DRadialCOM.C
    polyDensityMOB2DRadialCOMIO.C

\*---------------------------------------------------------------------------*/

#ifndef polyDensityMOB2DRadialCOM_H
#define polyDensityMOB2DRadialCOM_H

#include "polyField.H"
#include "vector.H"
#include "polyMesh.H"
#include "volFields.H"
#include "polyMoleculeCloud.H"
#include "distribution.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class polyDensityMOB2DRadialCOM Declaration
\*---------------------------------------------------------------------------*/

class polyDensityMOB2DRadialCOM
:
	public polyField
{

private:

    // Private data

        //- subDictionary containing the properties
        dictionary propsDict_;

        propertyField fields_;

        word fieldName_;

        word regionName_;
        label regionId_;
	
        // added by Jun Zhang for polymolecules	
	    List<label> molIds_;

        vector startPoint_;
        vector endPoint_;

        scalar h_;

        vector unitVectorX_;

        scalar radius_;

        label nBinsX_;
        label nBinsY_;

        scalar binWidthX_;
        scalar binWidthY_;

        Field<scalarField> mols_;
        vector oldCenterOfMass_;
        Field<scalarField> densityField_;

        scalarField magRadii_;
        scalarField binWidths_;
        scalarField volume_;

        scalar avVolume_;

        // for searching
        scalar minBinWidth_;
        labelField n_;
//         vector COM_;
        scalar nAvTimeSteps_;
        bool resetAtOutput_;
    // Private Member Functions

        label findBin(const scalar& r);


public:

    //- Runtime type information
    TypeName("polyDensityMOB2DRadialCOM");

    // Constructors

        //- Construct from components
        polyDensityMOB2DRadialCOM
        (
            Time& t,
            const polyMesh& mesh,
            polyMoleculeCloud& molCloud,
            const dictionary& dict
        );

    // Destructor

        ~polyDensityMOB2DRadialCOM();


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

//         void measureDuringForceComputation
//         (
//             polyMolecule* molReal,
//             polyReferredMolecule* molRef
//         );

        void measureDuringForceComputationSite
        (
            polyMolecule* molI,
            polyMolecule* molJ,
            label sI,
            label sJ
        );

//         void measureDuringForceComputationSite
//         (
//             polyMolecule* molReal,
//             polyReferredMolecule* molRef,
//             label sReal,
//             label sRef
//         );


        //- return field
//         const volVectorField& velocityField() const;
        const propertyField& fields() const;

        //- write field
        void writeField();

        void updateProperties(const dictionary&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
