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

#include "diffusivityModel.H"
#include "volFields.H"

namespace Foam
{
    defineTypeNameAndDebug(diffusivityModel, 0);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::diffusivityModel::diffusivityModel
(
    const word transportPropertiesDictName,
    const word thermoSpeciesDictName,
    const volScalarField& p,
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
    T_(T),
    species_(species)
    
{
    DijModels_.setSize(0.5*species.size()*(species.size()+1));
    Dij_.setSize(DijModels_.size());
           
    forAll(species, i)
    {
        for(label j=i; j < species.size(); j++)
        {
            label k = species.size()*i+j-0.5*i*(i+1);
            
            DijModels_.set
            (
                k,
                binaryDiffusivityModel::New
                (
                    species[i],
                    species[j],
                    dictThermo_,
                    dictTransport_,
                    p,
                    T
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

void Foam::diffusivityModel::update()
{    
    for(label i=0; i < species_.size(); i++)
    {
        for(label j=i; j < species_.size(); j++)
        {
            label k = species_.size()*i+j-0.5*i*(i+1);
            Dij_[k] = DijModels_[k].D();
        }
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
