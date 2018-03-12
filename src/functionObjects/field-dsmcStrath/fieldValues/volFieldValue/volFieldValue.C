/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  |
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

#include "volFieldValue.H"
#include "fvMesh.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
namespace fieldValues
{
    defineTypeNameAndDebug(volFieldValue, 0);
    addToRunTimeSelectionTable(fieldValue, volFieldValue, dictionary);
    addToRunTimeSelectionTable(functionObject, volFieldValue, dictionary);
}
}
}

template<>
const char*
Foam::NamedEnum
<
    Foam::functionObjects::fieldValues::volFieldValue::operationType,
    11
>::names[] =
{
    "none",
    "sum",
    "sumMag",
    "average",
    "weightedAverage",
    "volAverage",
    "weightedVolAverage",
    "volIntegrate",
    "min",
    "max",
    "CoV"
};

const Foam::NamedEnum
<
    Foam::functionObjects::fieldValues::volFieldValue::operationType,
    11
> Foam::functionObjects::fieldValues::volFieldValue::operationTypeNames_;


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::functionObjects::fieldValues::volFieldValue::initialise
(
    const dictionary& dict
)
{
    if (dict.readIfPresent("weightField", weightFieldName_))
    {
        Info<< "    weight field = " << weightFieldName_;
    }

    Info<< nl << endl;
}


void Foam::functionObjects::fieldValues::volFieldValue::writeFileHeader
(
    Ostream& os
) const
{
    volRegion::writeFileHeader(*this, os);

    writeCommented(os, "Time");

    forAll(fields_, fieldi)
    {
        os  << tab << operationTypeNames_[operation_]
            << "(" << fields_[fieldi] << ")";
    }

    os  << endl;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::fieldValues::volFieldValue::volFieldValue
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fieldValue(name, runTime, dict, typeName),
    volRegion(fieldValue::mesh_, dict),
    operation_(operationTypeNames_.read(dict.lookup("operation"))),
    weightFieldName_("none")
{
    read(dict);
    writeFileHeader(file());
}


Foam::functionObjects::fieldValues::volFieldValue::volFieldValue
(
    const word& name,
    const objectRegistry& obr,
    const dictionary& dict
)
:
    fieldValue(name, obr, dict, typeName),
    volRegion(fieldValue::mesh_, dict),
    operation_(operationTypeNames_.read(dict.lookup("operation"))),
    weightFieldName_("none")
{
    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::fieldValues::volFieldValue::~volFieldValue()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::fieldValues::volFieldValue::read
(
    const dictionary& dict
)
{
    fieldValue::read(dict);

    // No additional info to read
    initialise(dict);

    return true;
}


bool Foam::functionObjects::fieldValues::volFieldValue::write()
{
    fieldValue::write();

    if (Pstream::master())
    {
        writeTime(file());
    }

    // Construct weight field. Note: zero size means weight = 1
    scalarField weightField;
    if (weightFieldName_ != "none")
    {
        weightField =
            getFieldValues<scalar>
            (
                weightFieldName_,
                true
            );
    }

    forAll(fields_, i)
    {
        const word& fieldName = fields_[i];
        bool ok = false;

        ok = ok || writeValues<scalar>(fieldName, weightField);
        ok = ok || writeValues<vector>(fieldName, weightField);
        ok = ok || writeValues<sphericalTensor>(fieldName, weightField);
        ok = ok || writeValues<symmTensor>(fieldName, weightField);
        ok = ok || writeValues<tensor>(fieldName, weightField);

        if (!ok)
        {
            WarningInFunction
                << "Requested field " << fieldName
                << " not found in database and not processed"
                << endl;
        }
    }

    if (Pstream::master())
    {
        file()<< endl;
    }

    Log << endl;

    return true;
}


// ************************************************************************* //
