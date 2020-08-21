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

#include "VTRelaxationModel.H"
#include "volFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(VTRelaxationModel, 0);
    defineRunTimeSelectionTable(VTRelaxationModel, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::VTRelaxationModel::VTRelaxationModel
(
    const word& name1,
    const word& name2,
    const label& lname1,
    const label& lname2,
    const dictionary& dict1,
    const dictionary& dict2,
    const volScalarField& p,
    const volScalarField& T,
    const PtrList<volScalarField>& Tv,
    const PtrList<volScalarField>& nD
)
:
    dict1_(dict1),
    dict2_(dict2),
    name1_(name1),
    name2_(name2),
    lname1_(lname1),
    lname2_(lname2),
    p_(p),
    T_(T),
    Tv_(Tv),
    nD_(nD),
    VTFullCoeffsForm_
    (
        readBool
        (
            dict1_.subDict("thermalRelaxationModels").subDict("VT")
                .lookup("fullCoeffsForm")
        )
    ),
    VTOverwriteDefault_
    (
        readBool
        (
            dict1_.subDict("thermalRelaxationModels").subDict("VT")
                .lookup("overwriteDefault")
        )
    ),
    VTSpeciesDependent_
    (
        readBool
        (
            dict1_.subDict("thermalRelaxationModels").subDict("VT")
                .lookup("speciesDependent")
        )
    ),
    VTCollidingPartner_
    (
        readBool
        (
            dict1_.subDict("thermalRelaxationModels").subDict("VT")
                .lookup("collidingPair")
        )
    )
{}


// * * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::VTRelaxationModel> Foam::VTRelaxationModel::New
(
    const word& name1,
    const word& name2,
    const label& lname1,
    const label& lname2,
    const dictionary& dict1,
    const dictionary& dict2,
    const volScalarField& p,
    const volScalarField& T,
    const PtrList<volScalarField>& Tv,
    const PtrList<volScalarField>& nD
)
{
    word VTRelaxationModelTypeName
    (
        dict1.subDict("thermalRelaxationModels").subDict("VT").lookup("model")
    );

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(VTRelaxationModelTypeName);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "VTRelaxationModel::New(const volVectorField&, "
            "const surfaceScalarField&)"
        )   << "Unknown VTRelaxationModel type "
            << VTRelaxationModelTypeName << endl << endl
            << "Valid VTRelaxationModels are: " << endl
            << dictionaryConstructorTablePtr_->toc()
            << exit(FatalError);
    }

    return autoPtr<VTRelaxationModel>
        (cstrIter()(name1, name2, lname1, lname2, dict1, dict2, p, T, Tv, nD));
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
