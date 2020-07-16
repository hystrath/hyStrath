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
    delFromCylindricalShell

Description

SourceFiles
    delFromCylindricalShellI.H
    delFromCylindricalShell.C
    delFromCylindricalShellIO.C

\*---------------------------------------------------------------------------*/

#ifndef delFromCylindricalShell_H
#define delFromCylindricalShell_H

#include "molsToDeleteModel.H"
#include "vector.H"
#include "polyMesh.H"
#include "volFields.H"
#include "dsmcCloud.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class delFromCylindricalShell Declaration
\*---------------------------------------------------------------------------*/

class delFromCylindricalShell
:
	public molsToDeleteModel
{

private:

    // Private data

        //- subDictionary containing the properties
        dictionary propsDict_;

        vector startPoint_;
        vector endPoint_;
        vector unitVector_;

        scalar rOut_;
        scalar rIn_;

        labelList typeIds_;

    // Private Member Functions

        void findMolsToDel();

public:

    //- Runtime type information
    TypeName("delFromCylindricalShell");

    // Constructors

        //- Construct from components
        delFromCylindricalShell
		(
            dsmcCloud& cloud,
            const dictionary& dict
		);

    // Destructor

 		~delFromCylindricalShell();


    // Member Functions


};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
