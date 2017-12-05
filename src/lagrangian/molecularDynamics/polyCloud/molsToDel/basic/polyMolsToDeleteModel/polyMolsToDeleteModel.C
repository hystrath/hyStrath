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

#include "polyMolsToDeleteModel.H"
#include "IFstream.H"
#include "graph.H"
#include "polyMoleculeCloud.H"

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(polyMolsToDeleteModel, 0);
defineRunTimeSelectionTable(polyMolsToDeleteModel, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
polyMolsToDeleteModel::polyMolsToDeleteModel
(
    polyMoleculeCloud& molCloud,
    const dictionary& dict
)
:
    mesh_(refCast<const fvMesh>(molCloud.mesh())),
    molCloud_(molCloud)
{}

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

autoPtr<polyMolsToDeleteModel> polyMolsToDeleteModel::New
(
    polyMoleculeCloud& molCloud,
    const dictionary& dict
)
{
    word polyMolsToDeleteModelName
    (
        dict.lookup("model")
    );

    Info<< "Selecting model "
         << polyMolsToDeleteModelName << endl;

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(polyMolsToDeleteModelName);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalError
            << "polyMolsToDeleteModel::New(const dictionary&) : " << endl
            << "    unknown molsToDeleteModel type "
            << polyMolsToDeleteModelName
            << ", constructor not in hash table" << endl << endl
            << "    Valid types are :" << endl;
        Info<< dictionaryConstructorTablePtr_->toc() << abort(FatalError);
    }

    return autoPtr<polyMolsToDeleteModel>
	(
		cstrIter()(molCloud, dict)
	);
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

polyMolsToDeleteModel::~polyMolsToDeleteModel()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void polyMolsToDeleteModel::deleteMolFromMoleculeCloud
(
    polyMolecule& mol
)
{
    //- remove polyMolecule from cloud
    molCloud_.deleteParticle(mol);
}

} // End namespace Foam

// ************************************************************************* //
