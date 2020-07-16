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

#include "rdfModel.H"


namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(rdfModel, 0);

defineRunTimeSelectionTable(rdfModel, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
rdfModel::rdfModel
(
//     Time& t,
    const dictionary& dict
)
:
//     time_(t)
    dict_(dict)
{}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

autoPtr<rdfModel> rdfModel::New
(
//     Time& t,
    const dictionary& dict
)
{
    word rdfModelName
    (
        dict.lookup("rdfModel")
    );

    Info<< "Selecting rdf function model: "
         << rdfModelName << endl;

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(rdfModelName);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalError
            << "rdfModel::New(const dictionary&) : " << endl
            << "    unknown rdfModel type "
            << rdfModelName
            << ", constructor not in hash table" << endl << endl
            << "    Valid types are :" << endl;
        Info<< dictionaryConstructorTablePtr_->toc() << abort(FatalError);
    }

    return autoPtr<rdfModel>
	(
		cstrIter()(/*t,*/ dict)
	);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

rdfModel::~rdfModel()
{}

} // End namespace Foam

// ************************************************************************* //
