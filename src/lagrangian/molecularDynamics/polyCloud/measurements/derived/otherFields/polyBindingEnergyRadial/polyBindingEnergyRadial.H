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
    polyBindingEnergyRadial


Description
    a rectangular measurement box divided into nBinsX and nBinsY
    but rotated about an axis to measure rotationally symmetric averages



SourceFiles
    polyBindingEnergyRadialI.H
    polyBindingEnergyRadial.C
    polyBindingEnergyRadialIO.C

\*---------------------------------------------------------------------------*/

#ifndef polyBindingEnergyRadial_H
#define polyBindingEnergyRadial_H

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
                           Class polyBindingEnergyRadial Declaration
\*---------------------------------------------------------------------------*/

class polyBindingEnergyRadial
:
	public polyField
{

private:

    // Private data

        //- subDictionary containing the properties
        dictionary propsDict_;

        propertyField fields_;
        
        fileName pathName_;
        fileName nameFile_;
        
        word fieldName_;

        vector startPoint_; // only y component is important here
//         vector endPoint_;
        
        vector unitVectorX_;
        vector unitVectorY_; // height
        vector unitVectorZ_;        
        
        label nBinsX_;
        label nBinsY_; // height 

        scalar lengthX_; // radius
        scalar lengthY_; 
        
        scalar binWidthX_; // radius
        scalar binWidthY_; // height

        
        List<label> molIdsWall_;
        List<label> molIdsFluid_;

        vector centreOfMass_;

        scalarField binsX_;
        scalarField binsY_;

        List<scalarField> volumes_;
        
        List<scalarField> pairs_;
        List<scalarField> mols_;
        List<scalarField> energy_;        

        // fields
        List<scalarField> bindingEnergy_;
        List<scalarField> avNoOfPairs_;

        List<vector> molPositions_;
        List<bool> fluidMols_;
        List<scalar> nPairs_;
        List<scalar> DeltaE_;
        
        
//         scalar nAvTimeSteps_;
        
        bool resetAtOutput_;        
        
        
        
    // Private Member Functions

        bool readFromStorage();
        
        void writeToStorage();
        
        List<label> isPointWithinBin
        (
            const vector& rI
        );
        
        void centreOfMass();
        
public:

    //- Runtime type information
    TypeName("polyBindingEnergyRadial");

    // Constructors

        //- Construct from components
        polyBindingEnergyRadial
		(
            Time& t,
            const polyMesh& mesh,
            polyMoleculeCloud& molCloud,
            const dictionary& dict
		);

    // Destructor

 		~polyBindingEnergyRadial();


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

        const propertyField& fields() const;

        //- write field
        void writeField();
}; 		


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
