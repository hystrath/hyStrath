/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012-2016 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2015-2016 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

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

\*---------------------------------------------------------------------------*/

#include "fieldValueDelta.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
namespace fieldValues
{
    defineTypeNameAndDebug(fieldValueDelta, 0);
    addToRunTimeSelectionTable(functionObject, fieldValueDelta, dictionary);
}
}
}

template<>
const char* Foam::NamedEnum
<
    Foam::functionObjects::fieldValues::fieldValueDelta::operationType,
    5
>::names[] =
{
    "add",
    "subtract",
    "min",
    "max",
    "average"
};
const Foam::NamedEnum
<
    Foam::functionObjects::fieldValues::fieldValueDelta::operationType,
    5
> Foam::functionObjects::fieldValues::fieldValueDelta::operationTypeNames_;


// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

void Foam::functionObjects::fieldValues::fieldValueDelta::writeFileHeader
(
    Ostream& os
) const
{
    const wordList& fields1 = region1Ptr_->fields();
    const wordList& fields2 = region2Ptr_->fields();

    DynamicList<word> commonFields(fields1.size());
    forAll(fields1, fieldi)
    {
        label index = findIndex(fields2, fields1[fieldi]);
        if (index != -1)
        {
            commonFields.append(fields1[fieldi]);
        }
    }

    writeHeaderValue(os, "Source1", region1Ptr_->name());
    writeHeaderValue(os, "Source2", region2Ptr_->name());
    writeHeaderValue(os, "Operation", operationTypeNames_[operation_]);
    writeCommented(os, "Time");

    forAll(commonFields, i)
    {
        os  << tab << commonFields[i];
    }

    os  << endl;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::fieldValues::fieldValueDelta::fieldValueDelta
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    writeFile(obr_, name, typeName, dict),
    operation_(opSubtract),
    region1Ptr_(nullptr),
    region2Ptr_(nullptr)
{
    read(dict);
    writeFileHeader(file());
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::fieldValues::fieldValueDelta::~fieldValueDelta()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::fieldValues::fieldValueDelta::read
(
    const dictionary& dict
)
{
    fvMeshFunctionObject::read(dict);
    writeFile::read(dict);

    region1Ptr_.reset
    (
        fieldValue::New
        (
            name() + ".region1",
            obr_,
            dict.subDict("region1"),
            false
        ).ptr()
    );
    region2Ptr_.reset
    (
        fieldValue::New
        (
            name() + ".region2",
            obr_,
            dict.subDict("region2"),
            false
        ).ptr()
    );

    operation_ = operationTypeNames_.read(dict.lookup("operation"));

    return true;
}


bool Foam::functionObjects::fieldValues::fieldValueDelta::write()
{
    region1Ptr_->write();
    region2Ptr_->write();

    writeTime(file());

    Log << type() << " " << name() << " write:" << endl;

    const word& name1 = region1Ptr_->name();
    const word& name2 = region2Ptr_->name();

    const wordList entries1 = objectResultEntries(name1);
    const wordList entries2 = objectResultEntries(name2);

    if (entries1.size() != entries2.size())
    {
        FatalErrorInFunction
            << name() << ": objects must generate the same number of results"
            << nl
            << "    " << name1 << " objects: " << entries1 << nl
            << "    " << name2 << " objects: " << entries2 << nl
            << exit(FatalError);
    }

    forAll(entries1, i)
    {
        const word& entry1(entries1[i]);
        const word& entry2(entries2[i]);
        const word type1 = objectResultType(name1, entry1);
        const word type2 = objectResultType(name2, entry2);

        if (type1 != type2)
        {
            FatalErrorInFunction
                << name()
                << ": input values for operation must be of the same type"
                << nl
                << "    " << entry1 << ": " << type1 << nl
                << "    " << entry2 << ": " << type2 << nl
                << exit(FatalError);
        }

        bool found = false;

        applyOperation<scalar>(type1, name1, name2, entry1, entry2, found);
        applyOperation<vector>(type1, name1, name2, entry1, entry2, found);
        applyOperation<sphericalTensor>
            (type1, name1, name2, entry1, entry2, found);
        applyOperation<symmTensor>(type1, name1, name2, entry1, entry2, found);
        applyOperation<tensor>(type1, name1, name2, entry1, entry2, found);

        if (!found)
        {
            Log << "Operation between "
                << name1 << " with result " << entry1 << " and "
                << name2 << " with result " << entry2 << " not applied"
                << endl;
        }
    }

    Log << (entries1.empty() ? "    none" : "") << endl;

    file()<< endl;

    return true;
}


bool Foam::functionObjects::fieldValues::fieldValueDelta::execute()
{
    return true;
}


// ************************************************************************* //
