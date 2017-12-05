/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2005 OpenCFD Ltd.
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

Class
    gravityForce

Description

\*----------------------------------------------------------------------------*/

#include "gravityForce.H"
#include "graph.H"
#include "IFstream.H"
#include "OFstream.H"

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(gravityForce, 0);

defineRunTimeSelectionTable(gravityForce, dictionary);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

gravityForce::gravityForce
(
    Time& time,
    const dictionary& dict
)
:
    time_(time),
    spaceVarying_(false),
    timeVarying_(false)
{}

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

autoPtr<gravityForce> gravityForce::New
(
    Time& time,
    const dictionary& dict
)
{
    word gravityForceName
    (
        dict.lookup("model")
    );

    Info<< "Selecting gravity-force model "
         << gravityForceName << endl;

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(gravityForceName);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalError
            << "gravityForce::New(const dictionary&) : " << endl
            << "    unknown gravityForce type "
            << gravityForceName
            << ", constructor not in hash table" << endl << endl
            << "    Valid gravityForce types are :" << endl;
        Info<< dictionaryConstructorTablePtr_->toc() << abort(FatalError);
    }

    return autoPtr<gravityForce>
    (
        cstrIter()(time, dict)
    );
}
// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

gravityForce::~gravityForce()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool gravityForce::spaceVarying()
{
    return spaceVarying_;
}

bool gravityForce::timeVarying()
{
    return timeVarying_;
}

// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
