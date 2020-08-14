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

#include "VTModel.H"
#include "volFields.H"

namespace Foam
{
    defineTypeNameAndDebug(VTModel, 0);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::VTModel::VTModel
(
    const word& dict2T,
    const word& dictThermoPhy,
    const wordList& solvedVibEqSpecies,
    const wordList& species,
    const volScalarField& p,
    const volScalarField& Tt,
    const PtrList<volScalarField>& Tv,
    const PtrList<volScalarField>& nD
)
:
    dict2T_
    (
        IOobject
        (
            dict2T,
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
            dictThermoPhy,
            Tt.mesh().time().constant(),
            Tt.mesh(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),

    solvedVibEqSpecies_(solvedVibEqSpecies),
    species_(species),
    p_(p),
    T_(Tt),
    Tv_(Tv),
    nD_(nD)

{
    tauVTijModels_.setSize(species.size()*species.size());
    tauVTij_.setSize(tauVTijModels_.size());

    forAll(solvedVibEqSpecies_, i)
    {
        forAll(species_, j)
        {
            label k = solvedVibEqSpecies_.size()*i+j;

            tauVTijModels_.set
            (
                k,
                VTRelaxationModel::New
                (
                    solvedVibEqSpecies[i],
                    species[j],
                    i,
                    j,
                    dict2T_,
                    dictThermoPhy_,
                    p,
                    Tt,
                    Tv,
                    nD
                )
            );

            tauVTij_.set
            (
                k,
                new volScalarField
                (
                   tauVTijModels_[k].tauVT()
                )
            );
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::VTModel::update()
{
    for(label i=0; i < solvedVibEqSpecies_.size(); i++)
    {
        for(label j=0; j < species_.size(); j++)
        {
            label k = solvedVibEqSpecies_.size()*i+j;
            tauVTij_[k] = tauVTijModels_[k].tauVT();
        }
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
