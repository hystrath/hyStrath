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

#include "VVModel.H"
#include "volFields.H"

namespace Foam
{
    defineTypeNameAndDebug(VVModel, 0);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::VVModel::VVModel
(
    const word& dict1,
    const word& dict2,
    const wordList& species,
    const volScalarField& p,
    const volScalarField& Tt,
    const PtrList<volScalarField>& Tv,
    const PtrList<volScalarField>& nD
)
:
    dict1_
    (
        IOobject
        (
            dict1,
            Tt.mesh().time().constant(),
            Tt.mesh(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),

    dict2_
    (
        IOobject
        (
            dict2,
            Tt.mesh().time().constant(),
            Tt.mesh(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),

    species_(species),
    p_(p),
    T_(Tt),
    Tv_(Tv),
    nD_(nD)

{
    tauVVijModels_.setSize(species.size()*species.size());
    tauVVij_.setSize(tauVVijModels_.size());

    forAll(species, i)
    {
        forAll(species, j)
        {
            label k = species_.size()*i+j;

            tauVVijModels_.set
            (
                k,
                VVRelaxationModel::New
                (
                    species[i],
                    species[j],
                    i,
                    j,
                    dict1_,
                    dict2_,
                    p,
                    Tt,
                    Tv,
                    nD
                )
            );

            tauVVij_.set
            (
                k,
                new volScalarField
                (
                   tauVVijModels_[k].tauVV()
                )
            );
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::VVModel::update()
{
    for(label i=0; i < species_.size(); i++)
    {
        for(label j=0; j < species_.size(); j++)
        {
            label k = species_.size()*i+j;
            tauVVij_[k] = tauVVijModels_[k].tauVV();
        }
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
