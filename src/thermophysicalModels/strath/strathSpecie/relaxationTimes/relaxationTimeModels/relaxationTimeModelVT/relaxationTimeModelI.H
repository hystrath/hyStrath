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
Foam::relaxationTimeModel::tauVTij(label i, label j) const
{
    return tauVTijModel_().tauVT(i,j);
}


inline const Foam::volScalarField&
Foam::relaxationTimeModel::QVT(const label speciei) const
{
    return QVT_[speciei];
}


inline const Foam::volScalarField&
Foam::relaxationTimeModel::QVT(const word& specieName) const
{
    return QVT_[solvedVibEqSpecies()[specieName]];
}


inline Foam::volScalarField& Foam::relaxationTimeModel::QVT(const label speciei)
{
    return QVT_[speciei];
}


inline Foam::volScalarField&
Foam::relaxationTimeModel::QVT(const word& specieName)
{
    return QVT_[solvedVibEqSpecies()[specieName]];
}


/*inline const Foam::volScalarField& // TODO ABORTIVE WORK
Foam::relaxationTimeModel::QVT_mode(label speciei, label vibMode) const
{
    return QVTmode_[speciei][vibMode];
}*/


inline const Foam::speciesTable& Foam::relaxationTimeModel::species() const
{
    return thermo_.composition().species();
}


inline const Foam::speciesTable&
Foam::relaxationTimeModel::solvedVibEqSpecies() const
{
    return thermo_.composition().solvedVibEqSpecies();
}


// ************************************************************************* //
