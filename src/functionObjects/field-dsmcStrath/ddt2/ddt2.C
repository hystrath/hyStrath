/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016 OpenCFD Ltd.
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

#include "ddt2.H"

#include "volFields.H"
#include "dictionary.H"
#include "wordReListMatcher.H"
#include "steadyStateDdtScheme.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(ddt2, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        ddt2,
        dictionary
    );
}
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::functionObjects::ddt2::checkFormatName(const word& str)
{
    if (str.find("@@") == string::npos)
    {
        WarningInFunction
            << "Bad result naming "
            << "(no '@@' token found), deactivating."
            << nl << endl;

        return false;
    }
    else if (str == "@@")
    {
        WarningInFunction
            << "Bad result naming "
            << "(only a '@@' token found), deactivating."
            << nl
            << endl;

        return false;
    }
    else
    {
        return true;
    }
}


bool Foam::functionObjects::ddt2::accept(const word& fieldName) const
{
    // check input vs possible result names
    // to avoid circular calculations
    return !blacklist_.match(fieldName);
}


int Foam::functionObjects::ddt2::process(const word& fieldName)
{
    if (!accept(fieldName))
    {
        return -1;
    }

    int state = 0;

    apply<volScalarField>(fieldName, state);
    apply<volVectorField>(fieldName, state);

    return state;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::ddt2::ddt2
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    selectFields_(),
    resultName_(word::null),
    blacklist_(),
    results_(),
    mag_(dict.lookupOrDefault<Switch>("mag", false))
{
    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::ddt2::~ddt2()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::ddt2::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);

    if (word(mesh_.ddtScheme("default")) == "steadyState")
    {
        WarningInFunction
            << typeName << " function object not appropriate for steady-state"
            << endl;
        return false;
    }

    selectFields_ = wordReListMatcher::uniq
    (
        wordReList(dict.lookup("fields"))
    );
    Info<< type() << " fields: " << selectFields_ << nl;

    resultName_ = dict.lookupOrDefault<word>
    (
        "result",
        ( mag_ ? "mag(ddt(@@))" : "magSqr(ddt(@@))" )
    );

    if (checkFormatName(resultName_))
    {
        blacklist_.set
        (
            string::quotemeta<regExp>
            (
                resultName_
            ).replace("@@", "(.+)")
        );

        return true;
    }
    else
    {
        blacklist_.clear();
        return false;
    }
}


bool Foam::functionObjects::ddt2::execute()
{
    results_.clear();

    wordHashSet candidates = subsetStrings(selectFields_, mesh_.names());
    DynamicList<word> missing(selectFields_.size());
    DynamicList<word> ignored(selectFields_.size());

    // check exact matches first
    forAll(selectFields_, i)
    {
        const wordRe& select = selectFields_[i];
        if (!select.isPattern())
        {
            const word& fieldName = static_cast<const word&>(select);

            if (!candidates.erase(fieldName))
            {
                missing.append(fieldName);
            }
            else if (process(fieldName) < 1)
            {
                ignored.append(fieldName);
            }
        }
    }

    forAllConstIter(wordHashSet, candidates, iter)
    {
        process(iter.key());
    }

    if (missing.size())
    {
        WarningInFunction
            << "Missing field " << missing << endl;
    }
    if (ignored.size())
    {
        WarningInFunction
            << "Unprocessed field " << ignored << endl;
    }

    return true;
}


bool Foam::functionObjects::ddt2::write()
{
    if (results_.size())
    {
        Log << type() << ' ' << name() << " write:" << endl;
    }

    // Consistent output order
    const wordList outputList = results_.sortedToc();
    forAll(outputList, i)
    {
        const word& fieldName = outputList[i];

        if (foundObject<regIOobject>(fieldName))
        {
            const regIOobject& io = lookupObject<regIOobject>(fieldName);

            Log << "    " << fieldName << endl;

            io.write();
        }
    }

    return true;
}


// ************************************************************************* //
