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

Description

\*---------------------------------------------------------------------------*/

#include "molsToDelete.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

//- Constructor
molsToDelete::molsToDelete
(
    const polyMesh& mesh,
    dsmcCloud& cloud
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
	modelList_(molsToDeleteDict_.lookup("deleteMols")),
	delModels_(modelList_.size())
{
    if( delModels_.size() > 0 )
    {
        forAll(delModels_, dM)
        {
            const entry& molsToDeleteI = modelList_[dM];
            const dictionary& molsToDeleteIDict = molsToDeleteI.dict();

            Info << nl << "Deleting molecules from model #: " << dM << endl;

            delModels_[dM] = autoPtr<molsToDeleteModel>
            (
                molsToDeleteModel::New(cloud, molsToDeleteIDict)
            );
        }
    }
}




// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
