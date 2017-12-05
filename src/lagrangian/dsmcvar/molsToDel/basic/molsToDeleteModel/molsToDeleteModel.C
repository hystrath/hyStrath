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

#include "molsToDeleteModel.H"
#include "IFstream.H"
#include "graph.H"
#include "dsmcCloud.H"

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(molsToDeleteModel, 0);

defineRunTimeSelectionTable(molsToDeleteModel, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
molsToDeleteModel::molsToDeleteModel
(
    dsmcCloud& cloud,
    const dictionary& dict
)
:
    mesh_(cloud.mesh()),
    cloud_(cloud)
//     molsToDel_()
{
}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

autoPtr<molsToDeleteModel> molsToDeleteModel::New
(
    dsmcCloud& cloud,
    const dictionary& dict
)
{
    word molsToDeleteModelName
    (
        dict.lookup("molsToDeleteModel")
    );

    Info<< "Selecting molsToDeleteModel "
         << molsToDeleteModelName << endl;

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(molsToDeleteModelName);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalError
            << "molsToDeleteModel::New(const dictionary&) : " << endl
            << "    unknown molsToDeleteModel type "
            << molsToDeleteModelName
            << ", constructor not in hash table" << endl << endl
            << "    Valid types are :" << endl;
        Info<< dictionaryConstructorTablePtr_->toc() << abort(FatalError);
    }

    return autoPtr<molsToDeleteModel>
	(
		cstrIter()(cloud, dict)
	);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

molsToDeleteModel::~molsToDeleteModel()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// const labelList& molsToDeleteModel::molsToDel() const
// {
//     return molsToDel_;
// }


void molsToDeleteModel::deleteMolFromMoleculeCloud
(
    dsmcParcel& mol
)
{
    //- remove molecule from cloud
    cloud_.deleteParticle(mol);
}


} // End namespace Foam

// ************************************************************************* //
