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

\*---------------------------------------------------------------------------*/

#include "basicMfpModel.H"
#include "volFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(basicMfpModel, 0);
    defineRunTimeSelectionTable(basicMfpModel, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::basicMfpModel::basicMfpModel
(
    const word& name,
    const label& speciesIndex,
    const dictionary& dict,
    const dictionary& dictThermoPhy,
    const volScalarField& p,
    const volScalarField& Tt
)
:
    dict_(dict),
    dictThermoPhy_(dictThermoPhy),
    name_(name),
    speciesIndex_(speciesIndex),
    p_(p),
    T_(Tt)
{}


// * * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::basicMfpModel> Foam::basicMfpModel::New
(
    const word& name,
    const label& speciesIndex,
    const dictionary& dict,
    const dictionary& dictThermoPhy,
    const volScalarField& p,
    const volScalarField& Tt
)
{
    word basicMfpModelTypeName(dict.subDict("rarefiedParameters").lookup("mfpModel"));

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(basicMfpModelTypeName);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "basicMfpModel::New(const volVectorField&, "
            "const surfaceScalarField&)"
        )   << "Unknown basicMfpModel type "
            << basicMfpModelTypeName << endl << endl
            << "Valid  basicMfpModels are : " << endl
            << dictionaryConstructorTablePtr_->toc()
            << exit(FatalError);
    }

    return autoPtr<basicMfpModel>
        (cstrIter()(name, speciesIndex, dict, dictThermoPhy, p, Tt));
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
