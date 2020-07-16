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
    polyMolsToDelete

Description



\*---------------------------------------------------------------------------*/

#ifndef polyMolsToDelete_H
#define polyMolsToDelete_H

#include "List.H"
#include "IOdictionary.H"
#include "autoPtr.H"
#include "polyMesh.H"
#include "polyMolsToDeleteModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class polyMolsToDelete Declaration
\*---------------------------------------------------------------------------*/

class polyMolsToDelete
{
    // Private data

        //- The entire dictionary (containing multiple subDictionaries)
        IOdictionary polyMolsToDeleteDict_;

        //- state polyMolsToDelete
		PtrList<entry> modelList_;
        List< autoPtr<polyMolsToDeleteModel> > delModels_;

public:

    // Constructors

        //- Constructor
        polyMolsToDelete
        (
            const polyMesh& mesh,
            polyMoleculeCloud& molCloud
        );

    // Member Functions

        // Access
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#endif

// ************************************************************************* //
