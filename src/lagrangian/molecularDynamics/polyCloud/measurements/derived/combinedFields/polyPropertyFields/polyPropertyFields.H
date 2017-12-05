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
    polyPropertyFields

Description

    Geometric fields for standardised fields

    One specie or all species handled only.

SourceFiles
    polyPropertyFieldsI.H
    polyPropertyFields.C
    polyPropertyFieldsIO.C

\*---------------------------------------------------------------------------*/

#ifndef polyPropertyFields_H
#define polyPropertyFields_H

#include "polyField.H"
#include "vector.H"
#include "polyMesh.H"
#include "volFields.H"
#include "polyMoleculeCloud.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class polyPropertyFields Declaration
\*---------------------------------------------------------------------------*/

class polyPropertyFields
:
	public polyField
{

private:

    // Private data

        dictionary propsDict_;

        propertyField fields_;

        word fieldName_;

        volScalarField N_;
        volScalarField rhoN_;
        volScalarField rhoM_;
        volVectorField USAM_;
        volVectorField UCAM_;
        volScalarField T_;
        volScalarField p_;
        volTensorField stress_;

        scalarField mols_;
        scalarField mass_;
        vectorField mom_;
        vectorField velocityB_;
        scalarField kE_;
        scalarField dof_;
        scalarField angularKe_;
        tensorField kineticTensor_;
        tensorField virialTensor_;

        scalarField molsV_;
        scalarField massV_;
        vectorField momV_;
        vectorField velocity_;
        vectorField angularSpeed_;
        vectorField angularVelocity_;

        List<label> molIds_;

        scalar nAvTimeSteps_;
        bool resetAtOutput_;

    // Private Member Functions

public:

    //- Runtime type information
    TypeName("polyPropertyFields");

    // Constructors

        //- Construct from components
        polyPropertyFields
		(
            Time& t,
            const polyMesh& mesh,
            polyMoleculeCloud& molCloud,
            const dictionary& dict
		);

    // Destructor

 		~polyPropertyFields();

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
