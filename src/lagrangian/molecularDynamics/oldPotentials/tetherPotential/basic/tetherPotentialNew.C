/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

#include "tetherPotential.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::tetherPotential> Foam::tetherPotential::New
(
    const word& name,
    const reducedUnits& rU,
    const dictionary& tetherPotentialProperties
)
{
    word
        tetherPotentialTypeName
        (
            tetherPotentialProperties.lookup("tetherPotential")
        );

    Info<< nl << "Selecting tether potential "
        << tetherPotentialTypeName << " for "
        << name << endl;

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(tetherPotentialTypeName);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "tetherPotential::New()"
        )   << "Unknown tetherPotential type "
            << tetherPotentialTypeName << nl << nl
            << "Valid  tetherPotentials are: " << nl
            << dictionaryConstructorTablePtr_->toc()
            << exit(FatalError);
    }

    return autoPtr<tetherPotential>
        (cstrIter()(name, rU, tetherPotentialProperties));
}


// ************************************************************************* //
