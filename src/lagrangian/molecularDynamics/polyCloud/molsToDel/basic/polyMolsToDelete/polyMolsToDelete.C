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

Description

\*---------------------------------------------------------------------------*/

#include "polyMolsToDelete.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

//- Constructor
polyMolsToDelete::polyMolsToDelete
(
    const polyMesh& mesh,
    polyMoleculeCloud& molCloud
)
:
    polyMolsToDeleteDict_
    (
        IOobject
        (
            "molsToDeleteDict",
            mesh.time().system(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        )
    ),
	modelList_(polyMolsToDeleteDict_.lookup("poly")),
	delModels_(modelList_.size())
{
    if( delModels_.size() > 0 )
    {
        forAll(delModels_, dM)
        {
            const entry& polyMolsToDeleteI = modelList_[dM];
            const dictionary& polyMolsToDeleteIDict = polyMolsToDeleteI.dict();

            Info << nl << "Deleting polyMolecules from model #: " << dM << endl;

            delModels_[dM] = autoPtr<polyMolsToDeleteModel>
            (
                polyMolsToDeleteModel::New(molCloud, polyMolsToDeleteIDict)
            );
        }
    }
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
