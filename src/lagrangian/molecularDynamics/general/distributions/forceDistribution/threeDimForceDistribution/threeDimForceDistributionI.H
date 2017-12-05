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

inline const word& threeDimForceDistribution::name() const
{
    return name_;
}

inline const word& threeDimForceDistribution::filename() const
{
    return fileName_;
}

inline const vector& threeDimForceDistribution::startPoint() const
{
    return startPoint_;
}

inline const vector& threeDimForceDistribution::unitVectorX() const
{
    return unitVectorX_;
}

inline const vector& threeDimForceDistribution::unitVectorY() const
{
    return unitVectorY_;
}

inline const vector& threeDimForceDistribution::unitVectorZ() const
{
    return unitVectorZ_;
}


inline const scalar& threeDimForceDistribution::binWidthX() const
{
    return binWidthX_;
}

inline const scalar& threeDimForceDistribution::binWidthY() const
{
    return binWidthY_;
}

inline const scalar& threeDimForceDistribution::binWidthZ() const
{
    return binWidthZ_;
}


inline const label& threeDimForceDistribution::nBinsX() const
{
    return nBinsX_;
}

inline const label& threeDimForceDistribution::nBinsY() const
{
    return nBinsY_;
}

inline const label& threeDimForceDistribution::nBinsZ() const
{
    return nBinsZ_;
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
