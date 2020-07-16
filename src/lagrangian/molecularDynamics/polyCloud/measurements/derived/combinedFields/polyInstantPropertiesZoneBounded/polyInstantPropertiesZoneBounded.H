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
    polyInstantPropertiesZoneBounded


Description
    Author: Srinivasa

SourceFiles
    polyInstantPropertiesZoneBoundedI.H
    polyInstantPropertiesZoneBounded.C
    polyInstantPropertiesZoneBoundedIO.C

\*---------------------------------------------------------------------------*/

#ifndef polyInstantPropertiesZoneBounded_H
#define polyInstantPropertiesZoneBounded_H

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
                           Class polyInstantPropertiesZoneBounded Declaration
\*---------------------------------------------------------------------------*/

class polyInstantPropertiesZoneBounded
:
    public polyField
{

private:

    // Private data

        //- subDictionary containing the properties
        dictionary propsDict_;

        propertyField fields_;

        word fieldName_;

        List<boundedBox> boxes_;

        List<label> molIds_;

        DynamicList<vector> velocityField_;
//         DynamicList<vector> forceField_;

        // Private Member Functions

        void setBoundBoxes();



public:

	//- Runtime type information
	TypeName("polyInstantPropertiesZoneBounded");

	// Constructors

	//- Construct from components
	polyInstantPropertiesZoneBounded
	(
		Time& t,
		const polyMesh& mesh,
		polyMoleculeCloud& molCloud,
		const dictionary& dict
		);

	// Destructor

	~polyInstantPropertiesZoneBounded();


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

	//void updateProperties(const dictionary&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
