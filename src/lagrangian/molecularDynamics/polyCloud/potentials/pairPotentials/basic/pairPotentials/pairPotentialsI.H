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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


inline List< autoPtr<pairPotentialModel> >& pairPotentials::pairPots()
{
    return pairPotentials_;
}

inline const List< autoPtr<pairPotentialModel> >& pairPotentials::pairPots() const
{
    return pairPotentials_;
}

inline const List<word>& pairPotentials::pairNames() const
{
    return pairNames_;
}

inline const List<label>& pairPotentials::pairIds() const
{
    return pairIds_;
}

inline autoPtr<pairPotentialModel>& pairPotentials::electrostaticPot()
{
    return electrostaticPotential_;
}

inline const autoPtr<pairPotentialModel>& pairPotentials::electrostaticPot() const
{
    return electrostaticPotential_;
}



// inline autoPtr<pairPotentialModel>& pairPotentials::pairPot
// (
//     const label pairSiteA,
//     const label pairSiteB
// )
// {
//      label k = pairPotLinks_[pairSiteA][pairSiteB];
//
//      return pairPotentials_[k];
// }

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
