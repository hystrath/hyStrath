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

#include "polyField.H"



namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(polyField, 0);

defineRunTimeSelectionTable(polyField, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
polyField::polyField
(
    Time& t,
    const polyMesh& mesh,
    polyMoleculeCloud& molCloud,
    const dictionary& dict
)
:
    mesh_(refCast<const fvMesh>(mesh)),
    molCloud_(molCloud),
    time_(t),
    casePath_(),
    timePath_(),
    measureInterForces_(false),
    measureInterForcesSites_(false)
{}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

autoPtr<polyField> polyField::New
(
    Time& t,
    const polyMesh& mesh,
    polyMoleculeCloud& molCloud,
    const dictionary& dict
)
{
    word polyFieldName
    (
        dict.lookup("fieldModel")
    );

    Info<< "Selecting field: "
         << polyFieldName << endl;

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(polyFieldName);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalError
            << "polyField::New(const dictionary&) : " << endl
            << "    unknown polyField type "
            << polyFieldName
            << ", constructor not in hash table" << endl << endl
            << "    Valid types are :" << endl;
        Info<< dictionaryConstructorTablePtr_->toc() << abort(FatalError);
    }

    return autoPtr<polyField>
	(
		cstrIter()(t, mesh, molCloud, dict)
	);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

polyField::~polyField()
{}


const fileName& polyField::casePath() const
{
    return casePath_;
}

fileName& polyField::casePath()
{
    return casePath_;
}

const fileName& polyField::timePath() const
{
    return timePath_;
}

fileName& polyField::timePath()
{
    return timePath_;
}

const bool& polyField::measureInterForces() const
{
    return measureInterForces_;
}

bool& polyField::measureInterForces()
{
    return measureInterForces_;
}

const bool& polyField::measureInterForcesSites() const
{
    return measureInterForcesSites_;
}

bool& polyField::measureInterForcesSites()
{
    return measureInterForcesSites_;
}
} // End namespace Foam

// ************************************************************************* //
