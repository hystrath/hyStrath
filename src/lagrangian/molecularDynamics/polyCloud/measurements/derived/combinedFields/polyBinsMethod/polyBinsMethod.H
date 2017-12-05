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
    polyBinsMethod

Description
    The method of bins measures 1D distribution properties localised in a
    user-specified zone.

    Parallelised.

    Any binModel may be applied. So far implemented the following:

        normal, cylindrical, multi sub-zones, cuboid 

    Choice of properties to output also available.

SourceFiles
    polyBinsMethodI.H
    polyBinsMethod.C
    polyBinsMethodIO.C

\*---------------------------------------------------------------------------*/

#ifndef polyBinsMethod_H
#define polyBinsMethod_H

#include "polyField.H"
#include "vector.H"
#include "polyMesh.H"
#include "volFields.H"
#include "polyMoleculeCloud.H"
#include "distribution.H"
#include "binModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class polyBinsMethod Declaration
\*---------------------------------------------------------------------------*/

class polyBinsMethod
:
	public polyField
{

private:

    // Private data

        //- subDictionary containing the properties
        dictionary propsDict_;

        fileName pathName_;
        fileName nameFile_;

        propertyField fields_;

        autoPtr<binModel> binModel_;

        word regionName_;
        label regionId_;
        word fieldName_;

        List<label> molIds_;

        scalarField mols_;
        scalarField mass_;
        vectorField mom_;
        vectorField velocityB_;
        scalarField kE_;
        scalarField angularKeSum_;
        scalarField dof_;
        tensorField kineticTensor_;
        tensorField virialTensor_;

        scalarField molsV_;
        scalarField massV_;
        vectorField momV_;
        vectorField velocity_;
        vectorField angularSpeed_;
        vectorField angularVelocity_;

        // fields
        scalarField N_;
        scalarField rhoN_;
        scalarField rhoM_;
        vectorField USAM_;
        vectorField UCAM_;
        scalarField T_;
        tensorField stress_;
        scalarField p_;
        scalarField pVir_;
        scalarField pKin_;

        List<bool> outputField_;

        scalar nAvTimeSteps_;
        
//         bool averagingAcrossManyRuns_;
        bool resetAtOutput_;
        
        
        
    // Private Member Functions
        
        
            void writeToStorage();

            bool readFromStorage();

public:

    //- Runtime type information
    TypeName("polyBinsMethod");

    // Constructors

        //- Construct from components
        polyBinsMethod
		(
            Time& t,
            const polyMesh& mesh,
            polyMoleculeCloud& molCloud,
            const dictionary& dict
		);

    // Destructor

 		~polyBinsMethod();


    // Member Functions

        //- create initial field
        void createField();

        //- calculate field
        void calculateField();

        //- return field
        const propertyField& fields() const;

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
}; 		


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
