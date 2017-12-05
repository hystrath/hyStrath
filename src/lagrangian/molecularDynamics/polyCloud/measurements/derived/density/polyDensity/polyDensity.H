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
    polyDensity

Description
            instantaneous density 

SourceFiles
    polyDensity.C


\*---------------------------------------------------------------------------*/

#ifndef polyDensity_H
#define polyDensity_H

#include "polyField.H"
#include "vector.H"
#include "polyMesh.H"
#include "volFields.H"
#include "polyMoleculeCloud.H"
#include "boundedBox.H" 

// #include "distribution.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class polyDensity Declaration
\*---------------------------------------------------------------------------*/

class polyDensity
:
	public polyField
{

private:

    // Private data

        //- subDictionary containing the properties
        dictionary propsDict_;
    
        word fieldName_;

        propertyField fields_;

        boundedBox bb_;

        //- ids of molecules
        List<label> molIds_;
    
        scalar volume_;
        
        DynamicList<scalar> densities_;
        
        void setBoundBox
        (
            const dictionary& propsDict,
            boundedBox& bb,
            const word& name
        );

    // Private Member Functions


 public:

    //- Runtime type information
    TypeName("polyDensity");

    // Constructors

        //- Construct from components
        polyDensity
		(
            Time& t,
            const polyMesh& mesh,
            polyMoleculeCloud& molCloud,
            const dictionary& dict
		);

    // Destructor

 		~polyDensity();


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

        //- write field
        void writeField();

        const propertyField& fields() const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
