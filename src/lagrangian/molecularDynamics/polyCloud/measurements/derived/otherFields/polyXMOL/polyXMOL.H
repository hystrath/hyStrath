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
    polyXMOL


Description
    to write out a pdb file for creating movies

SourceFiles
    polyXMOLI.H
    polyXMOL.C
    polyXMOLIO.C

\*---------------------------------------------------------------------------*/

#ifndef polyXMOL_H
#define polyXMOL_H

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
                           Class polyXMOL Declaration
\*---------------------------------------------------------------------------*/

class polyXMOL
:
	public polyField
{

private:

    // Private data

        //- subDictionary containing the properties
        dictionary propsDict_;

        propertyField fields_;

        List<label> molIds_;

        List<word> excludeSites_;

        word fieldName_;

        label iteration_;

        // PDB in zone        
        bool zone_;
        word regionName_;
        label regionId_;

    // Private Member Functions

            void write();

            void writeInZone(List<labelField>& molIds, List<vectorField>& sites);

            void writeInMesh(List<labelField>& molIds, List<vectorField>& sites);

public:

    //- Runtime type information
    TypeName("polyXMOL");

    // Constructors

        //- Construct from components
        polyXMOL
		(
            Time& t,
            const polyMesh& mesh,
            polyMoleculeCloud& molCloud,
            const dictionary& dict
		);

    // Destructor

 		~polyXMOL();

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
