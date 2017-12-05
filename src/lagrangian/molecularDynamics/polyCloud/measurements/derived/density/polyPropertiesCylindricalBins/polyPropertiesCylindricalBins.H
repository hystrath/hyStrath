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
    polyPropertiesCylindricalBins

Description
    The method of bins: measuring 2D distribution of velocity.

    Parallelised.

    Only applicable for cuboid or cubic domains or domains with a constant cross-sectional area.

    Two methods for measuring velocity are used, 
    see Tysanner and Garcia J. Comp. Phys. 2004.

SourceFiles
    polyPropertiesCylindricalBinsI.H
    polyPropertiesCylindricalBins.C
    polyPropertiesCylindricalBinsIO.C

\*---------------------------------------------------------------------------*/

#ifndef polyPropertiesCylindricalBins_H
#define polyPropertiesCylindricalBins_H

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
                           Class polyPropertiesCylindricalBins Declaration
\*---------------------------------------------------------------------------*/

class polyPropertiesCylindricalBins
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
        
        // for sampling the velocity
        Field<scalarField> mass_;
        Field<scalarField> momX_;
        Field<scalarField> momY_;
        Field<scalarField> sqrMomX_;
        Field<scalarField> sqrMomY_;
        
        // for samling the temperature;
        Field<scalarField> kE_;
        Field<scalarField> dof_;
        
        vector oldCenterOfMass_;
        
        Field<scalarField> densityField_;
        
        // for output the velocity and temperature
        Field<scalarField> velocityFieldX_;
        Field<scalarField> velocityFieldY_;
        Field<scalarField> temperatureField_;

        scalarField magRadii_;
        scalarField binWidths_;
        scalarField volume_;

        scalar avVolume_;

        // for searching
        scalar minBinWidth_;
        labelField n_;
//         vector COM_;

    // Private Member Functions

        label findBin(const scalar& r);


public:

    //- Runtime type information
    TypeName("polyPropertiesCylindricalBins");

    // Constructors

        //- Construct from components
        polyPropertiesCylindricalBins
        (
            Time& t,
            const polyMesh& mesh,
            polyMoleculeCloud& molCloud,
            const dictionary& dict
        );

    // Destructor

        ~polyPropertiesCylindricalBins();


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
