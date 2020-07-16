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

#include "molsToDelete.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

//- Constructor
molsToDelete::molsToDelete
(
    const polyMesh& mesh,
    dsmcCloud& cloud,
    const word& deleteType
)
:
    molsToDeleteDict_
    (
        IOobject
        (
            "molsToDeleteDict",
            mesh.time().system(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE,
            false
        )
    ),
	  modelList_(molsToDeleteDict_.lookup("deleteMols" + deleteType)),
	  delModels_(modelList_.size())
{
    if (delModels_.size() > 0)
    {
        forAll(delModels_, dM)
        {
            const entry& molsToDeleteI = modelList_[dM];
            const dictionary& molsToDeleteIDict = molsToDeleteI.dict();

            Info<< nl << "Deleting particles from model " << deleteType
                << ": #" << dM << endl;

            delModels_[dM] = autoPtr<molsToDeleteModel>
            (
                molsToDeleteModel::New(cloud, molsToDeleteIDict)
            );
        }
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

molsToDelete::~molsToDelete()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void molsToDelete::update()
{
    forAll(delModels_, dM)
    {
        delModels_[dM]->update();
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
