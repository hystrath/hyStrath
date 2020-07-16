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
    threeDimBinModel

Description

\*----------------------------------------------------------------------------*/

#include "threeDimBinModel.H"
#include "graph.H"


namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(threeDimBinModel, 0);

defineRunTimeSelectionTable(threeDimBinModel, dictionary);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

threeDimBinModel::threeDimBinModel
(
    const polyMesh& mesh,
    const dictionary& dict
)
:
    mesh_(mesh)
{}

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

autoPtr<threeDimBinModel> threeDimBinModel::New
(
    const polyMesh& mesh,
    const dictionary& dict
)
{
    word threeDimBinModelName
    (
        dict.lookup("threeDimBinModel")
    );

    Info<< "Selecting threeDimBinModel "
         << threeDimBinModelName << endl;

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(threeDimBinModelName);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalError
            << "threeDimBinModel::New(const dictionary&) : " << endl
            << "    unknown threeDimBinModel type "
            << threeDimBinModelName
            << ", constructor not in hash table" << endl << endl
            << "    Valid threeDimBinModel types are :" << endl;
        Info<< dictionaryConstructorTablePtr_->toc() << abort(FatalError);
    }

    return autoPtr<threeDimBinModel>
    (
        cstrIter()(mesh, dict)
    );
}
// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

threeDimBinModel::~threeDimBinModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //



// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
