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

inline const Foam::volScalarField&
Foam::multiSpeciesTransportModel::Dij(label i, label j) const
{
    return DijModel_().D(i,j);
}


inline const Foam::volScalarField&
Foam::multiSpeciesTransportModel::rhoD(label i) const
{
    return D_[i];
}


inline const Foam::volScalarField&
Foam::multiSpeciesTransportModel::Di(label i) const
{
    return rhoD(i)/thermo_.rho();
}


inline const Foam::surfaceScalarField&
Foam::multiSpeciesTransportModel::speciesMassFlux(const label i) const
{
    return spMassFlux_[i];
}


inline const Foam::volVectorField&
Foam::multiSpeciesTransportModel::J(const label i) const
{
    return JnonCorrected_[i];
}


inline const Foam::volVectorField&
Foam::multiSpeciesTransportModel::JGradp(const label i) const
{
    return JGradp_[i];
}


inline const Foam::volVectorField&
Foam::multiSpeciesTransportModel::JGradT(const label i) const
{
    return JGradT_[i];
}


inline const Foam::volVectorField
Foam::multiSpeciesTransportModel::Jcorr(const label i) const
{
    return Jcorrected(i);
}


inline const Foam::volVectorField&
Foam::multiSpeciesTransportModel::sumDiffusiveFluxes() const
{
    return sumDiffusionFluxes_;
}


inline const Foam::volVectorField&
Foam::multiSpeciesTransportModel::diffusionVelocity(const label i) const
{
    return Jcorrected(i)/thermo_.rho();
}


inline const Foam::speciesTable&
Foam::multiSpeciesTransportModel::species() const
{
    return thermo_.composition().species();
}

// ************************************************************************* //
