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

inline const word& twoDimForceDistribution::name() const
{
    return name_;
}

inline const word& twoDimForceDistribution::filename() const
{
    return fileName_;
}

inline const vector& twoDimForceDistribution::startPoint() const
{
    return startPoint_;
}

inline const vector& twoDimForceDistribution::unitVectorX() const
{
    return unitVectorX_;
}

inline const vector& twoDimForceDistribution::unitVectorY() const
{
    return unitVectorY_;
}

inline const scalar& twoDimForceDistribution::binWidthX() const
{
    return binWidthX_;
}

inline const scalar& twoDimForceDistribution::binWidthY() const
{
    return binWidthY_;
}

inline const label& twoDimForceDistribution::nBinsX() const
{
    return nBinsX_;
}

inline const label& twoDimForceDistribution::nBinsY() const
{
    return nBinsY_;
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
// inline List< Pair<scalar> > radialDistribution::radialDistributionFunction () const
// {
//     return radialDistribution_;
// }

// inline List< scalar > radialDistribution::radialDistributionFunction () const
// {
//     return radialDistribution_;
// }
// 
// 
// inline List< scalar > radialDistribution::radialPosition() const
// {
//     return radialPosition_;
// }

// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
