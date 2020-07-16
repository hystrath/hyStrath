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

#include "dsmcField.H"
#include "graph.H"


namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(dsmcField, 0);

defineRunTimeSelectionTable(dsmcField, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
dsmcField::dsmcField
(
    Time& t,
    const polyMesh& mesh,
    dsmcCloud& cloud,
    const dictionary& dict
)
:
    mesh_(refCast<const fvMesh>(mesh)),
    cloud_(cloud),
    timeDict_(dict.subDict("timeProperties")),
    time_(t, timeDict_),
    casePath_(),
    timePath_()
{
}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

autoPtr<dsmcField> dsmcField::New
(
    Time& t,
    const polyMesh& mesh,
    dsmcCloud& cloud,
    const dictionary& dict
)
{
    word dsmcFieldName
    (
        dict.lookup("fieldModel")
    );

    Info<< "Selecting field: "
         << dsmcFieldName << endl;

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(dsmcFieldName);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalError
            << "dsmcField::New(const dictionary&) : " << endl
            << "    unknown dsmcField type "
            << dsmcFieldName
            << ", constructor not in hash table" << endl << endl
            << "    Valid types are :" << endl;
        Info<< dictionaryConstructorTablePtr_->toc() << abort(FatalError);
    }

    return autoPtr<dsmcField>
	(
		cstrIter()(t, mesh, cloud, dict)
	);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

dsmcField::~dsmcField()
{}

void dsmcField::updateTime()
{
    time_++;
}


void dsmcField::updateBasicFieldProperties
(
    const dictionary& newDict
)
{
    timeDict_ = newDict.subDict("timeProperties");

    if (timeDict_.found("resetAtOutput"))
    {
        time_.resetFieldsAtOutput() = Switch(timeDict_.lookup("resetAtOutput"));
    }
}

const fileName& dsmcField::casePath() const
{
    return casePath_;
}

fileName& dsmcField::casePath()
{
    return casePath_;
}


const fileName& dsmcField::timePath() const
{
    return timePath_;
}

fileName& dsmcField::timePath()
{
    return timePath_;
}



} // End namespace Foam

// ************************************************************************* //
