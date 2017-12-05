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
    polyForceTwoMolecules


Description
    Measures directly the force in a two-molecule simulation

SourceFiles
    polyForceTwoMoleculesI.H
    polyForceTwoMolecules.C
    polyForceTwoMoleculesIO.C

\*---------------------------------------------------------------------------*/

#ifndef polyForceTwoMolecules_H
#define polyForceTwoMolecules_H

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
                           Class polyForceTwoMolecules Declaration
\*---------------------------------------------------------------------------*/

class polyForceTwoMolecules
:
	public polyField
{

private:

    // Private data

        //- subDictionary containing the properties
        dictionary propsDict_;

        propertyField fields_;

        word fieldName_;

        // poly site id
        label molSiteIdA_;
        label molSiteIdB_;

        vector molPointA_;
        vector molPointB_;

        vector molPositionA_;
        vector molPositionB_;

        label cellA_;
        label cellB_;

        label molTrackingNoA_;
        label molTrackingNoB_;
        
        scalar writeInterval_;
        
        scalar rMax_;
        label nBins_;
        scalar binWidth_;
        scalarField r_;
        scalarField f_;

        tensorField virA_;
        tensorField virB_;

        scalarField mols_;

        label timeIndex_;
        scalarField separation_;

    // Private Member Functions

        void setMolIds();

public:

    //- Runtime type information
    TypeName("polyForceTwoMolecules");

    // Constructors

        //- Construct from components
        polyForceTwoMolecules
		(
            Time& t,
            const polyMesh& mesh,
            polyMoleculeCloud& molCloud,
            const dictionary& dict
		);

    // Destructor

 		~polyForceTwoMolecules();

    // Member Functions

        //- create initial field
        void createField();

        //- calculate field
        void calculateField();

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
