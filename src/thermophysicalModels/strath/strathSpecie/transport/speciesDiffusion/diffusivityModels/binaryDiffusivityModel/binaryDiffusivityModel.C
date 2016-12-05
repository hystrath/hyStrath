/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
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

\*---------------------------------------------------------------------------*/

#include "binaryDiffusivityModel.H"
#include "volFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(binaryDiffusivityModel, 0);
    defineRunTimeSelectionTable(binaryDiffusivityModel, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::binaryDiffusivityModel::binaryDiffusivityModel
(
    const word& name1,
    const word& name2,
    const dictionary& dictThermo,
    const dictionary& dictTransport,
    const volScalarField& p,
    const volScalarField& T
)
:
    name1_(name1),    
    name2_(name2),
    dictThermo_(dictThermo),
    dictTransport_(dictTransport),    
    p_(p),   
    T_(T)
{}


// * * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::binaryDiffusivityModel> Foam::binaryDiffusivityModel::New
(
    const word& name1,
    const word& name2,
    const dictionary& dictThermo,
    const dictionary& dictTransport,
    const volScalarField& p,
    const volScalarField& T
)
{
    word binaryDiffusivityModelTypeName(dictTransport.subDict("transportModels").lookup("binaryDiffusivityModel"));

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(binaryDiffusivityModelTypeName);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "DiffusivityModel::New(const volVectorField&, "
            "const surfaceScalarField&)"
        )   << "Unknown binaryDiffusivityModel type "
            << binaryDiffusivityModelTypeName << endl << endl
            << "Valid  binaryDiffusivityModels are : " << endl
            << dictionaryConstructorTablePtr_->toc()
            << exit(FatalError);
    }

    return autoPtr<binaryDiffusivityModel>
        (cstrIter()(name1, name2, dictThermo, dictTransport, p, T));
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
