/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
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

#include "eVModel.H"
#include "volFields.H"

namespace Foam
{
    defineTypeNameAndDebug(eVModel, 0);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::eVModel::eVModel
(
    const word& dict2T,
    const word& dictThermoPhy,
    const wordList& species,
    const volScalarField& p,
    const volScalarField& Tv
)
:
    dict2T_
    (
        IOobject
        (
            dict2T,
            Tv.mesh().time().constant(),
            Tv.mesh(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    
    dictThermoPhy_
    (
        IOobject
        (
            dictThermoPhy,
            Tv.mesh().time().constant(),
            Tv.mesh(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),

    species_(species), 
    p_(p),    
    Tv_(Tv)
{
    taueViModels_.setSize(species.size());
    taueVi_.setSize(taueViModels_.size());
    
    forAll(species, speciei)
    {
        taueViModels_.set
        (
            speciei,
            eVRelaxationModel::New
            (
                species[speciei],
                speciei,
                dict2T_,
                dictThermoPhy_,
                p,
                Tv
            )
        );

        taueVi_.set
        (
            speciei,
            new volScalarField
            (
               taueViModels_[speciei].taueV()
            )
        );
    } 
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::eVModel::update()
{    
    forAll(species_, speciei)
    {
        taueVi_[speciei] = taueViModels_[speciei].taueV();
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
