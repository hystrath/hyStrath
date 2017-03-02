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

inline const Foam::volScalarField&
Foam::relaxationTimeModeleV::taueVi(label i) const
{
    return taueViModel_().taueV(i);
}


inline const Foam::volScalarField&
Foam::relaxationTimeModeleV::QeV(label speciei) const
{
    return QeV_[speciei];
}


inline const Foam::speciesTable&
Foam::relaxationTimeModeleV::species() const
{
    return thermo_.composition().species();
}


inline const Foam::speciesTable&
Foam::relaxationTimeModeleV::solvedVibEqSpecies() const
{
    return thermo_.composition().solvedVibEqSpecies();
}

// ************************************************************************* //
