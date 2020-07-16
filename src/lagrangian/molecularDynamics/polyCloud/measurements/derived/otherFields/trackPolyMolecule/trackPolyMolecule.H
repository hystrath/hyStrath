/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016-2020 hyStrath
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of hyStrath, a derivative work of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class

Description
        Tracks an individual molecule through the mesh. Specify the approximate
        start point of the molecule and its id. The code will loop and find the
        closest molecule to the specified start point and store its "trackingNumber".
        At every averaging time instance of the simulation run, the code will
        loop over all molecules of the domain until it finds that molecule with
        the target trackingNumber and store its current position to a list.

        The code is parallelised.

        More work can be done to make this efficient. (avoid looping over all
        molecules in the domain)

SourceFiles
    trackPolyMoleculeI.H
    trackPolyMolecule.C
    trackPolyMoleculeIO.C

\*---------------------------------------------------------------------------*/

#ifndef trackPolyMolecule_H
#define trackPolyMolecule_H

#include "polyField.H"
#include "vector.H"
#include "polyMesh.H"
#include "volFields.H"
#include "polyMoleculeCloud.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class trackPolyMolecule Declaration
\*---------------------------------------------------------------------------*/

class trackPolyMolecule
:
	public polyField
{

private:

    // Private data

        //- subDictionary containing the properties
        dictionary propsDict_;

        propertyField fields_;

        word fieldName_;

        label molId_;

        vector startPoint_;

        label trackingNumber_;

        vectorField trackingOrientation1_;
        vectorField trackingOrientation2_;
        vectorField trackingOrientation3_;
        tensorField trackingQ_;
        vectorField trackingPi_;
	    vectorField trackingLinMom_;
        vectorField trackingPositions_;

        Field<vectorField> trackingSitePositions_;

        scalarField phi_;
        scalarField theta_;

        labelField cells_;

        bool outputParaFoamFile_;

        bool outputVMDFile_;

    // Private Member Functions


public:

    //- Runtime type information
    TypeName("trackPolyMolecule");

    // Constructors

        //- Construct from components
        trackPolyMolecule
		(
            Time& t,
            const polyMesh& mesh,
            polyMoleculeCloud& molCloud,
            const dictionary& dict
		);

    // Destructor

 		~trackPolyMolecule();

    // Member Functions

        //- create initial field
        void createField();

        //- calculate field
        void calculateField();

        //- write field
        void writeField();

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
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
