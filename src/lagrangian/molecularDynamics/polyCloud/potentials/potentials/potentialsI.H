/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

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

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

namespace Foam
{
    

inline const scalar& potentials::potentialEnergyLimit() const
{
    return potentialEnergyLimit_;
}
/*
inline Foam::label Foam::potential::nPairPotentials() const
{
    return pairPotentials_.size();
}*/


inline const labelList& potentials::removalOrder() const
{
    return removalOrder_;
}

/*
inline const Foam::pairPotentialList& Foam::potential::pairPotentials() const
{
    return pairPotentials_;
}*/
} // End namespace Foam

// ************************************************************************* //
