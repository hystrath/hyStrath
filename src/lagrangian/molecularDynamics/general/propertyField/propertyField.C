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

Description

\*---------------------------------------------------------------------------*/

#include "propertyField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

propertyField::propertyField
(
    Time& t,
    const polyMesh& mesh,
    const word& fieldName
)
:
    mesh_(refCast<const fvMesh>(mesh)),
    time_(t),
    fieldName_(fieldName),
    sField_
    (
        IOobject
        (
            fieldName_+"_scalar",
            time_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("0.0", dimless, 0.0)
    ),
    vField_
    (
        IOobject
        (
            fieldName_+"_vector",
            time_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedVector("zero", dimless, vector::zero)
    ),
    tField_
    (
        IOobject
        (
            fieldName_+"_tensor",
            time_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedTensor("zero", dimless, tensor::zero)
    ),
    s_(0.0),
    v_(vector::zero),
    t_(tensor::zero)
{}




propertyField::~propertyField()
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
