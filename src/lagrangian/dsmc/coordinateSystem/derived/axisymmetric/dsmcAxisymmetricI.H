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

Class

Description

\*---------------------------------------------------------------------------*/

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline scalar dsmcAxisymmetric::radialExtent() const
{
    return radialExtent_;
}


inline scalar dsmcAxisymmetric::maxRWF() const
{
    return maxRWF_;
}


inline const volScalarField& dsmcAxisymmetric::nParticles() const
{
    return timeStepModel_().nParticles()*RWF_;
}


inline scalar dsmcAxisymmetric::nParticles(const label celli) const
{
    return timeStepModel_().nParticles(celli)*RWF(celli);
}


inline scalar dsmcAxisymmetric::nParticles
(
    const label patchi,
    const label facei
) const
{
    return timeStepModel_().nParticles(patchi, facei)*pRWF(patchi, facei);
}


inline scalar dsmcAxisymmetric::RWF(const label celli) const
{
    return RWF_[celli];
}


inline scalar dsmcAxisymmetric::pRWF
(
    const label patchi,
    const label facei
) const
{
    return RWF_.boundaryField()[patchi][facei];
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
