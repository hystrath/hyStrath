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
    polyHydroxylSiO2SurfaceBounded

Description

SourceFiles
    polyHydroxylSiO2SurfaceBoundedI.H
    polyHydroxylSiO2SurfaceBounded.C
    polyHydroxylSiO2SurfaceBoundedIO.C

\*---------------------------------------------------------------------------*/

#ifndef polyHydroxylSiO2SurfaceBounded_H
#define polyHydroxylSiO2SurfaceBounded_H

#include "polyMolsToDeleteModel.H"
#include "vector.H"
#include "polyMesh.H"
#include "volFields.H"
#include "polyMoleculeCloud.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class polyHydroxylSiO2SurfaceBounded Declaration
\*---------------------------------------------------------------------------*/

class polyHydroxylSiO2SurfaceBounded
:
	public polyMolsToDeleteModel
{

private:

    // Private data

        //- subDictionary containing the properties
        dictionary propsDict_;

        
        // bond radius
        scalar b_;
        
        label SiId_;
        label OId_;
        label HId_;       
        
        boundBox box_;
        
        vector hDirection_;
        
        scalar temperature_;
        vector bulkVelocity_;
        
        scalar radius_;
        vector startPoint_;
        vector endPoint_; 
        scalar deltaR_;
        vector unitVector_;        
        
        DynamicList<label> siliconAtoms_;
        DynamicList<label> oxygenAtoms_; 
        
        List<label> siliconToOxygenBonds_;
        List<label> oxygenToSiliconBonds_;
        
    // Private Member Functions
        void readInCylinder();
        
        void readInBoundBox();
        
        void punchOutMolecules();
        
        void checkBoundBox
        (
            boundBox& b,
            const vector& startPoint,
            const vector& endPoint
        );
        
        void check();
      
        bool testForPair
        (
            const label& idI,
            const label& idJ,
            label& tnI,
            label& tnJ
        );
        
        void checkForMoleculesInBox();
        void checkForMoleculesInCylinder();
        
        void deleteMoleculesInBox();
        void deleteMoleculesInCylinder();
        
//         void addMoleculesInBox();
        
        void addMoleculesInBoxAndCylinder();
        
//         void insertMolecule
//         (
//             const point& position,
//             label& cell,
//             const label& id,
//             const bool& tethered,
//             const bool& frozen,
//             const scalar& temperature,
//             const vector& bulkVelocity
//         );
        
        vector equipartitionLinearVelocity
        (
            scalar temperature,
            scalar mass
        );

        vector equipartitionAngularMomentum
        (
            scalar temperature,
            label id
        );

public:

    //- Runtime type information
    TypeName("polyHydroxylSiO2SurfaceBounded");

    // Constructors

        //- Construct from components
        polyHydroxylSiO2SurfaceBounded
		(
            polyMoleculeCloud& molCloud,
            const dictionary& dict
		);

    // Destructor

 		~polyHydroxylSiO2SurfaceBounded();


    // Member Functions

 
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
