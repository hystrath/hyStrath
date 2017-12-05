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
    polyFrictionForce


Description
    output every write interval folder for this to work 

SourceFiles
    polyFrictionForceI.H
    polyFrictionForce.C
    polyFrictionForceIO.C

\*---------------------------------------------------------------------------*/

#ifndef polyFrictionForce_H
#define polyFrictionForce_H

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
                           Class polyFrictionForce Declaration
\*---------------------------------------------------------------------------*/

class polyFrictionForce
:
	public polyField
{

private:

    // Private data

        //- subDictionary containing the properties
        dictionary propsDict_;

        propertyField fields_;
        
//         fileName pathName_;
//         fileName nameFile_;
        
        word fieldName_;

        List<label> molIdsFluid_;
        List<label> molIdsWall_;
        
        vector startPoint_;
        vector endPoint_;
        
        // choose this small enough ~0.5 nm
        vector endPointS_;
        
        label nMolsMin_;
        
        scalar dr_;
        scalar radius_;
        label nBins_;
        scalar binWidth_;
        
        scalar h_;
        scalar hS_;
        
        vector n_;
        
        DynamicList<scalar> radii_;
        DynamicList<scalar> radii_;
        
        
        
        
        
        
//         List<scalar> volume_;
        
        void refreshBins();
        
        void setCentreOfMass();
        
        void setRadius();
        
        
    // Private Member Functions

//         bool readFromStorage();
        
//         void writeToStorage();
        
        List<label> isPointWithinBin
        (
            const vector& rI
        );
        
        void centreOfMass();
        
public:

    //- Runtime type information
    TypeName("polyFrictionForce");

    // Constructors

        //- Construct from components
        polyFrictionForce
		(
            Time& t,
            const polyMesh& mesh,
            polyMoleculeCloud& molCloud,
            const dictionary& dict
		);

    // Destructor

 		~polyFrictionForce();


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
