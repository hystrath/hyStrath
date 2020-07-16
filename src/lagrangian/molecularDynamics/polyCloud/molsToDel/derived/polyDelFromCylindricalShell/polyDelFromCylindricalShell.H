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
    polyDelFromCylindricalShell

Description

SourceFiles
    polyDelFromCylindricalShellI.H
    polyDelFromCylindricalShell.C
    polyDelFromCylindricalShellIO.C

\*---------------------------------------------------------------------------*/

#ifndef polyDelFromCylindricalShell_H
#define polyDelFromCylindricalShell_H

#include "polyMolsToDeleteModel.H"
#include "vector.H"
#include "polyMesh.H"
#include "volFields.H"
#include "polyMoleculeCloud.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class polyDelFromCylindricalShell Declaration
\*---------------------------------------------------------------------------*/

class polyDelFromCylindricalShell
:
	public polyMolsToDeleteModel
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

        List<label> molIds_;

    // Private Member Functions

        void findMolsToDel();

public:

    //- Runtime type information
    TypeName("polyDelFromCylindricalShell");

    // Constructors

        //- Construct from components
        polyDelFromCylindricalShell
		(
            polyMoleculeCloud& molCloud,
            const dictionary& dict
		);

    // Destructor

 		~polyDelFromCylindricalShell();


    // Member Functions


};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
