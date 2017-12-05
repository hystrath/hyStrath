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
    distribution

Description

\*---------------------------------------------------------------------------*/

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Destructors * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline const label& polyFadeII::nInserted() const
{
    return nInserted_;
}

inline const label& polyFadeII::nDeleted() const
{
    return nDeleted_;
}

inline const DynamicList<label>& polyFadeII::insertingList() const
{
    return trackingNumbersIns_;
}            

inline const DynamicList<label>& polyFadeII::deletingList() const
{
    return trackingNumbersDel_;
}

inline const List< DynamicList<vector> > & polyFadeII::overlappingMols() const
{
    return overlappingMolsPositions_;
}

inline bool polyFadeII::isNotFading(const label& tN)
{
    if
    (
        (findIndex(trackingNumbersIns_, tN) == -1)
        &&
        (findIndex(trackingNumbersDel_, tN) == -1)
    )
    {
        return true;
    }
    else
    {
        return false;
    }
}

inline const scalar& polyFadeII::temperature() const
{
    return temperature_;
}

inline scalar& polyFadeII::temperature()
{
    return temperature_;
}

inline const vector& polyFadeII::velocity() const
{
    return velocity_;
}

inline vector& polyFadeII::velocity()
{
    return velocity_;
}

// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
