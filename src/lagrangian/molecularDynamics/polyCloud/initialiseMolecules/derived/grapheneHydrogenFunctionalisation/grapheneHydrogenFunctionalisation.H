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
    grapheneHydrogenFunctionalisation

Description


SourceFiles
    grapheneHydrogenFunctionalisationI.H
    grapheneHydrogenFunctionalisation.C
    grapheneHydrogenFunctionalisationIO.C


\*---------------------------------------------------------------------------*/

#ifndef grapheneHydrogenFunctionalisation_H
#define grapheneHydrogenFunctionalisation_H

#include "polyConfiguration.H"
#include "vector.H"
#include "polyMesh.H"
#include "volFields.H"
#include "polyMoleculeCloud.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class grapheneHydrogenFunctionalisation Declaration
\*---------------------------------------------------------------------------*/

class grapheneHydrogenFunctionalisation
:
	public polyConfiguration
{

private:

    // Private data

        //- subDictionary containing the properties
//         dictionary propsDict_;
//         fileName timePath_;

        word option_;

        scalar dr_;

        label molIdC_;
        label molIdCN_;
        label molIdH_;

    // Private Member Functions


        void insertMolecule
        (
            const point& position,
            const label& cell,
            const label& id,
            const label& tetFace,
            const label& tetPt,
            const bool& tethered,
            const bool& frozen,
            const scalar& phi,
            const scalar& theta,
            const scalar& psi,
//             const scalar& temperature,
            const vector& velocity
        );

        void fixedPropertiesFromDict();



public:

    //- Runtime type information
    TypeName("grapheneHydrogenFunctionalisation");

    // Constructors

        //- Construct from components
        grapheneHydrogenFunctionalisation
		(
            polyMoleculeCloud& molCloud,
            const dictionary& dict
//             const word& name
		);

    // Destructor

 		~grapheneHydrogenFunctionalisation();


    // Member Functions

        void setInitialConfiguration();

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
