/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2005 OpenCFD Ltd.
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

Class

Description

\*---------------------------------------------------------------------------*/

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline scalar dsmcTimeStepModel::deltaTValueOrg() const
{
    return mesh_.time().deltaTValue();
}


inline const scalarField& dsmcTimeStepModel::deltaTField() const
{
    tmp<volScalarField> tdeltaT
    (
        new volScalarField
        (
            IOobject
            (
                "deltaT",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar("deltaT", dimTime, deltaTValueOrg())
        )
    );

    return tdeltaT.ref().primitiveField();
}


inline scalar dsmcTimeStepModel::deltaTValue(const label celli) const
{
    return deltaTValueOrg();
}


inline const volScalarField& dsmcTimeStepModel::nParticles() const
{
    return nParticles_;
}


inline scalar dsmcTimeStepModel::nParticles(const label celli) const
{
    return nParticles_[celli];
}


inline scalar dsmcTimeStepModel::nParticles
(
    const label patchi,
    const label facei
) const
{
    return nParticles_.boundaryField()[patchi][facei];
}


inline scalar dsmcTimeStepModel::nParticles(const label celli)
{
    return nParticles_[celli];
}


inline scalar dsmcTimeStepModel::nParticles
(
    const label patchi,
    const label facei
)
{
    return nParticles_.boundaryField()[patchi][facei];
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
