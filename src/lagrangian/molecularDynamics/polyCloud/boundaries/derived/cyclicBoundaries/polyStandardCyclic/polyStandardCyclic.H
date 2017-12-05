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
    polyStandardCyclic

Description


SourceFiles
    polyStandardCyclic.C


\*---------------------------------------------------------------------------*/

#ifndef polyStandardCyclic_H
#define polyStandardCyclic_H

#include "polyCyclicBoundary.H"
#include "vector.H"
#include "polyMesh.H"
#include "volFields.H"
#include "polyMoleculeCloud.H"
#include "distribution.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class polyStandardCyclic Declaration
\*---------------------------------------------------------------------------*/

class polyStandardCyclic
:
	public polyCyclicBoundary
{

private:

    // Private data
        dictionary propsDict_;



public:

    //- Runtime type information
    TypeName("polyStandardCyclic");

    // Constructors

        //- Construct from components
        polyStandardCyclic
		(
            Time& t,
            const polyMesh& mesh,
            polyMoleculeCloud& molCloud,
            const dictionary& dict
		);

    // Destructor

 		~polyStandardCyclic();


    // Member Functions

        void calculateProperties();

        //- create an initial configuration
        void initialConfiguration();

        //- control the molecules
        void controlMol
        (
            polyMolecule& mol,
            polyMolecule::trackingData& td,
            const label& patchi
        );

//         void correctAfterParallelTransfer
//         (
//             polyMolecule& mol,
//             const label& patchi,
//             polyMolecule::trackData& td
//         );

        void controlAfterMove();

        //- output data
        void output
        (
            const fileName& fixedPathName,
            const fileName& timePath
        );

        //- update properties from a modified dictionary
        void updateProperties(const dictionary&);




        // Access





};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
