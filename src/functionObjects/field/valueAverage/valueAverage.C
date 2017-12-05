/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015 OpenFOAM Foundation
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

#include "valueAverage.H"
#include "addToRunTimeSelectionTable.H"
#include "Time.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(valueAverage, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        valueAverage,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * //

void Foam::functionObjects::valueAverage::writeFileHeader(Ostream& os) const
{
    writeHeader(os, "Value averages");
    writeCommented(os, "Time");
    forAll(fieldNames_, fieldi)
    {
        writeTabbed(os, fieldNames_[fieldi]);
    }
    os  << endl;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::valueAverage::valueAverage
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    regionFunctionObject(name, runTime, dict),
    writeFile(obr_, name, typeName, dict),
    functionObjectName_("unknown-functionObject"),
    fieldNames_(),
    window_(-1),
    totalTime_(),
    resetOnRestart_(false)
{
    read(dict);

    if (resetOnRestart_)
    {
        forAll(fieldNames_, fieldi)
        {
            const word& fieldName = fieldNames_[fieldi];

            if (dict.found(fieldName))
            {
                const dictionary& valueDict = dict.subDict(fieldName);
                totalTime_[fieldi] = readScalar(valueDict.lookup("totalTime"));
            }
        }
    }

    writeFileHeader(file());
}


// * * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * //

Foam::functionObjects::valueAverage::~valueAverage()
{}


// * * * * * * * * * * * * * * Public Member Functions * * * * * * * * * * * //

bool Foam::functionObjects::valueAverage::read(const dictionary& dict)
{
    regionFunctionObject::read(dict);
    writeFile::read(dict);

    // Make certain that the values are consistent with the defaults:
    resetOnRestart_ = false;

    dict.lookup("functionObject") >> functionObjectName_;
    dict.lookup("fields") >> fieldNames_;
    window_ = dict.lookupOrDefault<scalar>("window", -1);

    totalTime_.setSize(fieldNames_.size());
    forAll(totalTime_, i)
    {
        totalTime_[i] = time_.deltaTValue();
    }

    dict.readIfPresent("resetOnRestart", resetOnRestart_);

    return true;
}


bool Foam::functionObjects::valueAverage::execute()
{
    scalar dt = obr_.time().deltaTValue();

    Log << type() << ": " << name() << " averages:" << nl;

    file() << time_.timeName();

    DynamicList<label> unprocessedFields(fieldNames_.size());

    forAll(fieldNames_, fieldi)
    {
        const word& fieldName(fieldNames_[fieldi]);
        const word meanName(fieldName + "Mean");

        scalar Dt = totalTime_[fieldi];
        scalar alpha = (Dt - dt)/Dt;
        scalar beta = dt/Dt;

        if (window_ > 0)
        {
            if (Dt - dt >= window_)
            {
                alpha = (window_ - dt)/window_;
                beta = dt/window_;
            }
        }

        bool processed = false;
        calc<scalar>(fieldName, meanName, alpha, beta, processed);
        calc<vector>(fieldName, meanName, alpha, beta, processed);
        calc<sphericalTensor>(fieldName, meanName, alpha, beta, processed);
        calc<symmTensor>(fieldName, meanName, alpha, beta, processed);
        calc<tensor>(fieldName, meanName, alpha, beta, processed);

        if (!processed)
        {
            unprocessedFields.append(fieldi);

            if (writeToFile())
            {
                file() << tab << "n/a";
            }
        }

        totalTime_[fieldi] += dt;
    }

    file()<< endl;

    if (unprocessedFields.size())
    {
        WarningInFunction
            << "From function object: " << functionObjectName_ << nl
            << "Unprocessed fields:" << nl;

        forAll(unprocessedFields, i)
        {
            label fieldi = unprocessedFields[i];
            Log << "        " << fieldNames_[fieldi] << nl;
        }
        Log << endl;
    }

    Log << endl;

    return true;
}


bool Foam::functionObjects::valueAverage::write()
{
    return true;
}


// ************************************************************************* //
