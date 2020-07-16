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
    wallMotion

Description

\*----------------------------------------------------------------------------*/

#include "wallMotion.H"
#include "graph.H"
#include "IFstream.H"
#include "OFstream.H"

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(wallMotion, 0);

defineRunTimeSelectionTable(wallMotion, dictionary);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

wallMotion::wallMotion
(
    Time& time,
    const dictionary& dict
)
:
    time_(time)
{}

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

autoPtr<wallMotion> wallMotion::New
(
    Time& time,
    const dictionary& dict
)
{
    word wallMotionName
    (
        dict.lookup("wallMotion")
    );

    Info<< "Selecting wallMotion "
         << wallMotionName << endl;

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(wallMotionName);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalError
            << "wallMotion::New(const dictionary&) : " << endl
            << "    unknown wallMotion type "
            << wallMotionName
            << ", constructor not in hash table" << endl << endl
            << "    Valid wallMotion types are :" << endl;
        Info<< dictionaryConstructorTablePtr_->toc() << abort(FatalError);
    }

    return autoPtr<wallMotion>
    (
        cstrIter()(time, dict)
    );
}
// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

wallMotion::~wallMotion()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //



// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
