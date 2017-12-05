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
