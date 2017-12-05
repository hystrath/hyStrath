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
    dsmcFreeStreamInflowFieldPatch

Description

    Fully deletion boundary patch: atomisticMolecules get deleted instantly upon
    collision with a face on the boundary. (aka vacuum boundary condition)

    Selective deletion of one specie is possible. Other species get
    specularly reflected.

SourceFiles
    dsmcFreeStreamInflowFieldPatch.C


\*---------------------------------------------------------------------------*/

#ifndef dsmcFreeStreamInflowFieldPatch_H
#define dsmcFreeStreamInflowFieldPatch_H

#include "dsmcGeneralBoundary.H"
#include "vector.H"
#include "polyMesh.H"
#include "volFields.H"
#include "dsmcCloud.H"
#include "distribution.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class dsmcFreeStreamInflowFieldPatch Declaration
\*---------------------------------------------------------------------------*/

class dsmcFreeStreamInflowFieldPatch
:
	public dsmcGeneralBoundary
{

private:

    // Private data
        dictionary propsDict_;

    labelList typeIds_;
    
    //vector so different translational, rotational
    //vibrational temperatures can be set
    vectorField inletTemperatures_; 
    
    vectorField inletVelocities_;

    // list of (one value) mass densities corresponding to the type Ids
    List <scalarField> numberDensities_;

    // list of species 
    List <scalarField> accumulatedParcelsToInsert_;

    //- boundary temperatures
    volVectorField boundaryT_;

    //- boundary velocity
    volVectorField boundaryU_;
    
    //- list of inlet densities (one entry for each species)
    List<autoPtr<volScalarField> > boundaryNumberDensity_;

    void setProperties();

public:

    //- Runtime type information
    TypeName("dsmcFreeStreamInflowFieldPatch");

    // Constructors

        //- Construct from components
        dsmcFreeStreamInflowFieldPatch
        (
            Time& t,
            const polyMesh& mesh,
            dsmcCloud& cloud,
            const dictionary& dict
        );

    // Destructor

        ~dsmcFreeStreamInflowFieldPatch();


    // Member Functions

        void initialConfiguration();

        //- calculate any required properties
        void calculateProperties();

        //- control the mols
        void controlParcelsBeforeMove();
        void controlParcelsBeforeCollisions();
        void controlParcelsAfterCollisions();

        //- output data
        void output
        (
            const fileName& fixedPathName,
            const fileName& timePath
        );

        //- update properties from a modified dictionary
        void updateProperties(const dictionary&);
        
        //reset fields if mesh has changed
        void setNewBoundaryFields();



        // Access

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
