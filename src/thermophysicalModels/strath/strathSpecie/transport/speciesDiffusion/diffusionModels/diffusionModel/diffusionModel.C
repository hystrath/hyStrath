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

#include "diffusionModel.H"
#include "volFields.H"

namespace Foam
{
    defineTypeNameAndDebug(diffusionModel, 0);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::diffusionModel::diffusionModel
(
    const word transportPropertiesDictName,
    const word thermoSpeciesDictName,
    const volScalarField& p,
    const volScalarField& pe,
    const volScalarField& T,
    const wordList& species
)
:
    dictTransport_
    (
        IOobject
        (
            transportPropertiesDictName,
            T.mesh().time().constant(),
            T.mesh(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),

    dictThermo_
    (
        IOobject
        (
            thermoSpeciesDictName,
            T.mesh().time().constant(),
            T.mesh(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),

    p_(p),
    pe_(pe),
    T_(T),
    species_(species)

{
    DijModels_.setSize(0.5*species.size()*(species.size()+1));
    Dij_.setSize(DijModels_.size());

    forAll(species, i)
    {
        for(label j=i; j < species.size(); j++)
        {
            const label k = species.size()*i+j-0.5*i*(i+1);

            DijModels_.set
            (
                k,
                binaryDiffusionModel::New
                (
                    species[i],
                    species[j],
                    dictThermo_,
                    dictTransport_,
                    p_,
                    pe_,
                    T_
                )
            );

            Dij_.set
            (
                k,
                new volScalarField
                (
                   DijModels_[k].D()
                )
            );
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::diffusionModel::update()
{
    forAll(species_, i)
    {
        for(label j=i; j < species_.size(); j++)
        {
            const label k = species_.size()*i+j-0.5*i*(i+1);
//            Info<< "k = " << i << j << endl;
            Dij_[k] = DijModels_[k].D();
//            Info<< "Dij_[k] = " << Dij_[k].internalField()[0] << endl;
        }
    }
    
//    FatalErrorIn
//    (
//        "DiffusionModel my error"
//    )   << exit(FatalError);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
