/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2009-2010 OpenCFD Ltd.
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

#include "ElectronModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


namespace Foam
{
    defineTypeNameAndDebug(ElectronModel, 0);

    defineRunTimeSelectionTable(ElectronModel, dictionary);
};



Foam::ElectronModel::ElectronModel(pdCloud& owner)
:
    dict_(dictionary::null),
    cloud_(owner)
{

}


//
Foam::ElectronModel::ElectronModel
(
    const dictionary& dict,
    pdCloud& owner
//     const word& type
)
:
    dict_(dict),
    cloud_(owner)
{
}

Foam::autoPtr<Foam::ElectronModel> Foam::ElectronModel::New
(
    const dictionary& dict,
    pdCloud& owner
)
{
    const word modelType(dict.lookup("ElectronModel"));

    Info<< "Selecting ElectronModel " << modelType << endl;

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(modelType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "ElectronModel::New"
            "(const dictionary&, CloudType&)"
        )
            << "Unknown ElectronModel type "
            << modelType << nl << nl
            << "Valid ElectronModel types are:" << nl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<ElectronModel>
    (
        cstrIter()(dict, owner)
    );
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //


Foam::ElectronModel::~ElectronModel()
{
    //Info << "ElectronModel Destructor" << endl;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::dictionary& Foam::ElectronModel::dict() const
{
    return dict_;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// #include "ElectronModelNew.C"

// ************************************************************************* //
