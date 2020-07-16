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
    distribution

Description

\*---------------------------------------------------------------------------*/

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Destructors * * * * * * * * * * * * * * * //

inline const word& forceDistribution::name() const
{
    return name_;
}

inline const word& forceDistribution::filename() const
{
    return fileName_;
}

inline const label& forceDistribution::noOfBins() const
{
    return noOfBins_;
}

inline const scalar& forceDistribution::binWidth() const
{
    return binWidth_;
}

inline const vector& forceDistribution::startPoint() const
{
    return startPoint_;
}

inline const vector& forceDistribution::endPoint() const
{
    return endPoint_;
}

inline const vector& forceDistribution::unitVector() const
{
    return unitVector_;
}

inline const vectorField& forceDistribution::radii() const
{
    return radii_;
}

inline const scalarField& forceDistribution::magRadii() const
{
    return magRadii_;
}

inline const vectorField& forceDistribution::forces() const
{
    return forces_;
}

inline vectorField& forceDistribution::forces()
{
    return forces_;
}

inline const scalarField& forceDistribution::energies() const
{
    return energies_;
}

inline scalarField& forceDistribution::energies()
{
    return energies_;
}


inline const scalarField& forceDistribution::magForces() const
{
//     vector rES = endPoint_ - startPoint_;
//
//     scalarField magForces = forces_ & (rES/mag(rES));

    return magForces_;
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
