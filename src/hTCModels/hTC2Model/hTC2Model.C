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

#include "hTC2Model.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(hTC2Model, 0);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::hTC2Model::hTC2Model
(
    const word& modelType,
    const fvMesh& mesh
)
:
    IOdictionary
    (
        IOobject
        (
            "hTCProperties",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    ),
    turbulencePtr_(),
    mesh_(mesh),
    active_(lookupOrDefault<Switch>("active", true)),
    coeffs_(subDict(modelType + "Coeffs")),
    modelType_(modelType)
{}


// * * * * * * * * * * * * * * * * Destructors * * * * * * * * * * * * * * * //

Foam::hTC2Model::~hTC2Model()
{
    if (turbulencePtr_)
    {
        turbulencePtr_ = 0;
    }
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //


bool Foam::hTC2Model::read()
{
    if (regIOobject::read())
    {
        this->lookup("active") >> active_;
        coeffs_ = subDict(modelType_ + "Coeffs");
        return true;
    }
    else
    {
        return false;
    }
}


Foam::tmp<Foam::volScalarField> Foam::hTC2Model::Sh() const
{
    return tmp<Foam::volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "Sh",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar("zero", dimEnergy/dimVolume/dimTime, 0.0)
        )
    );
}


Foam::tmp<Foam::volScalarField> Foam::hTC2Model::Scv() const
{
    return tmp<Foam::volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "Scv",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar("zero", dimEnergy/dimVolume/dimTime, 0.0)
        )
    );
}


Foam::tmp<Foam::volScalarField> Foam::hTC2Model::Scv(const label i) const
{
    return tmp<Foam::volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "Scv_" + word(i),
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar("zero", dimEnergy/dimVolume/dimTime, 0.0)
        )
    );
}


Foam::tmp<Foam::volScalarField> Foam::hTC2Model::Seiir() const
{
    return tmp<Foam::volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "Seiir",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar("zero", dimEnergy/dimVolume/dimTime, 0.0)
        )
    );
}


Foam::tmp<Foam::volScalarField> Foam::hTC2Model::Seiir(const label i) const
{
    return tmp<Foam::volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "Seiir_" + word(i),
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar("zero", dimEnergy/dimVolume/dimTime, 0.0)
        )
    );
}


// ************************************************************************* //
