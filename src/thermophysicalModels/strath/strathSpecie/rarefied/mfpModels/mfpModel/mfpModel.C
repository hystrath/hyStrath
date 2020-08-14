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

#include "mfpModel.H"
#include "volFields.H"

namespace Foam
{
    defineTypeNameAndDebug(mfpModel, 0);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::mfpModel::mfpModel
(
    const word transportPropertiesDictName,
    const word thermoSpeciesDictName,
    const wordList& species,
    const volScalarField& p,
    const volScalarField& Tt
)
:
    dict_
    (
        IOobject
        (
            transportPropertiesDictName,
            Tt.mesh().time().constant(),
            Tt.mesh(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),

    dictThermoPhy_
    (
        IOobject
        (
            thermoSpeciesDictName,
            Tt.mesh().time().constant(),
            Tt.mesh(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),

    species_(species),
    p_(p),
    T_(Tt)

{
    mfpModels_.setSize(species.size());
    mfp_.setSize(mfpModels_.size());

    forAll(species, i)
    {
        mfpModels_.set
        (
            i,
            basicMfpModel::New
            (
                species[i],
                i,
                dict_,
                dictThermoPhy_,
                p,
                Tt
            )
        );

        mfp_.set
        (
            i,
            new volScalarField
            (
               mfpModels_[i].mfp()
            )
        );
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::mfpModel::update()
{
    for(label i=0; i<species_.size(); i++)
    {
        mfp_[i] = mfpModels_[i].mfp();
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
