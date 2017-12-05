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
    poly2DBinsSites

Description


SourceFiles
    poly2DBinsSitesI.H
    poly2DBinsSites.C
    poly2DBinsSitesIO.C

\*---------------------------------------------------------------------------*/

#ifndef poly2DBinsSites_H
#define poly2DBinsSites_H

#include "polyField.H"
#include "vector.H"
#include "polyMesh.H"
#include "volFields.H"
#include "polyMoleculeCloud.H"
#include "distribution.H"
#include "twoDimBinModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class poly2DBinsSites Declaration
\*---------------------------------------------------------------------------*/

class poly2DBinsSites
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

        autoPtr<twoDimBinModel> binModel_;

        word regionName_;
        label regionId_;
        word fieldName_;

        List<label> molSiteIds_;

        List<scalarField> atoms_;


        // fields
        // av number of molecules
        List<scalarField> N_;


        scalar nAvTimeSteps_;
        
        bool resetAtOutput_;
        
        
        
    // Private Member Functions
        
        
            void writeToStorage();

            bool readFromStorage();

public:

    //- Runtime type information
    TypeName("poly2DBinsSites");

    // Constructors

        //- Construct from components
        poly2DBinsSites
		(
            Time& t,
            const polyMesh& mesh,
            polyMoleculeCloud& molCloud,
            const dictionary& dict
		);

    // Destructor

 		~poly2DBinsSites();


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
