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

#include "advancedSpecie.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

inline advancedSpecie::advancedSpecie
(
    const word& name,
    const scalar nMoles,
    const scalar molWeight,
    const label particleType,
    const label particleCharge,
    const scalar diameter,
    const scalar omega,
    const DynamicList<scalar> vibrationalList,
    const scalar dissociationPotential,
    const label noVibrationalTemp,
    const label noElectronicLevels,
    const scalar iHat
)
:
    name_(name),
    nMoles_(nMoles),
    molWeight_(molWeight),
    particleType_(particleType),
    particleCharge_(particleCharge),
    diameter_(diameter),
    omega_(omega),
    vibrationalList_(vibrationalList),
    dissociationPotential_(dissociationPotential),
    noVibrationalTemp_(noVibrationalTemp),
    noElectronicLevels_(noElectronicLevels),
    iHat_(iHat)
{}


inline advancedSpecie::advancedSpecie
(
    const scalar nMoles,
    const scalar molWeight,
    const label particleType,
    const label particleCharge,
    const scalar diameter,
    const scalar omega,
    const DynamicList<scalar> vibrationalList,
    const scalar dissociationPotential,
    const label noVibrationalTemp,
    const label noElectronicLevels,
    const scalar iHat
)
:
    nMoles_(nMoles),
    molWeight_(molWeight),
    particleType_(particleType),
    particleCharge_(particleCharge),
    diameter_(diameter),
    omega_(omega),
    vibrationalList_(vibrationalList),
    dissociationPotential_(dissociationPotential),
    noVibrationalTemp_(noVibrationalTemp),
    noElectronicLevels_(noElectronicLevels),
    iHat_(iHat)
{}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline advancedSpecie::advancedSpecie(const advancedSpecie& as)
:
    name_(as.name_),
    nMoles_(as.nMoles_),
    molWeight_(as.molWeight_),
    particleType_(as.particleType_),
    particleCharge_(as.particleCharge_),
    diameter_(as.diameter_),
    omega_(as.omega_),
    vibrationalList_(as.vibrationalList_),
    dissociationPotential_(as.dissociationPotential_),
    noVibrationalTemp_(as.noVibrationalTemp_),
    noElectronicLevels_(as.noElectronicLevels_),
    iHat_(as.iHat_)
{}


inline advancedSpecie::advancedSpecie(const word& name, const advancedSpecie& as)
:
    name_(name),
    nMoles_(as.nMoles_),
    molWeight_(as.molWeight_),
    particleType_(as.particleType_),
    particleCharge_(as.particleCharge_),
    diameter_(as.diameter_),
    omega_(as.omega_),
    vibrationalList_(as.vibrationalList_),
    dissociationPotential_(as.dissociationPotential_),
    noVibrationalTemp_(as.noVibrationalTemp_),
    noElectronicLevels_(as.noElectronicLevels_),
    iHat_(as.iHat_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline const word& advancedSpecie::nameSpecies() const
{
    return name_;
}


inline scalar advancedSpecie::W() const
{
    return molWeight_;
}


inline scalar advancedSpecie::nMoles() const
{
    return nMoles_;
}


inline scalar advancedSpecie::R() const
{
    return RR/molWeight_;
}


inline label advancedSpecie::particleType() const
{
    return particleType_;
}


inline label advancedSpecie::particleCharge() const
{
    return particleCharge_;
}


inline scalar advancedSpecie::diameter() const
{
    return diameter_;
}


inline scalar advancedSpecie::omega() const
{
    return omega_;
}


inline DynamicList<scalar> advancedSpecie::vibrationalList() const
{
    return vibrationalList_;
}


inline scalar advancedSpecie::dissociationPotential() const
{
    return dissociationPotential_;
}


inline label advancedSpecie::noVibrationalTemp() const
{
    return noVibrationalTemp_;
}


inline label advancedSpecie::noElectronicLevels() const
{
    return noElectronicLevels_;
}


inline scalar advancedSpecie::iHat() const
{
    return iHat_;
}

// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

inline void advancedSpecie::operator=(const advancedSpecie& as)
{
    nMoles_ = as.nMoles_;
    molWeight_ = as.molWeight_;
    particleType_ = as.particleType_;
    particleCharge_ = as.particleCharge_;
    diameter_ = as.diameter_;
    omega_ = as.omega_;
    vibrationalList_ = as.vibrationalList_;
    dissociationPotential_ = as.dissociationPotential_;
    noVibrationalTemp_ = as.noVibrationalTemp_;
    noElectronicLevels_ = as.noElectronicLevels_;
    iHat_ = as.iHat_;
}


inline void advancedSpecie::operator+=(const advancedSpecie& as)
{
    scalar sumNmoles = max(nMoles_ + as.nMoles_, SMALL);

    molWeight_ =
        nMoles_/sumNmoles*molWeight_
      + as.nMoles_/sumNmoles*as.molWeight_;

    nMoles_ = sumNmoles;
}


inline void advancedSpecie::operator-=(const advancedSpecie& as)
{
    scalar diffnMoles = nMoles_ - as.nMoles_;
    if (mag(diffnMoles) < SMALL)
    {
        diffnMoles = SMALL;
    }

    molWeight_ =
        nMoles_/diffnMoles*molWeight_
      - as.nMoles_/diffnMoles*as.molWeight_; 

    nMoles_ = diffnMoles;
}


inline void advancedSpecie::operator*=(const scalar s)
{
    nMoles_ *= s;
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

inline advancedSpecie operator+(const advancedSpecie& as1, const advancedSpecie& as2)
{
    scalar sumNmoles = max(as1.nMoles_ + as2.nMoles_, SMALL);

    return advancedSpecie
    (
        sumNmoles,
        as1.nMoles_/sumNmoles*as1.molWeight_
      + as2.nMoles_/sumNmoles*as2.molWeight_,
        as1.particleType_,
        as1.particleCharge_,
        as1.nMoles_/sumNmoles*as1.diameter_
      + as2.nMoles_/sumNmoles*as2.diameter_,
        as1.nMoles_/sumNmoles*as1.omega_
      + as2.nMoles_/sumNmoles*as2.omega_,
        as1.vibrationalList_,
        as1.nMoles_/sumNmoles*as1.dissociationPotential_
      + as2.nMoles_/sumNmoles*as2.dissociationPotential_,
        as1.noVibrationalTemp_,
        as1.noElectronicLevels_,
        as1.nMoles_/sumNmoles*as1.iHat_
      + as2.nMoles_/sumNmoles*as2.iHat_
    );
}


inline advancedSpecie operator-(const advancedSpecie& as1, const advancedSpecie& as2)
{
    scalar diffNmoles = as1.nMoles_ - as2.nMoles_;
    if (mag(diffNmoles) < SMALL)
    {
        diffNmoles = SMALL;
    }

    return advancedSpecie
    (
        diffNmoles,
        as1.nMoles_/diffNmoles*as1.molWeight_
      - as2.nMoles_/diffNmoles*as2.molWeight_,
        as1.particleType_,
        as1.particleCharge_,
        as1.nMoles_/diffNmoles*as1.diameter_
      - as2.nMoles_/diffNmoles*as2.diameter_,
        as1.nMoles_/diffNmoles*as1.omega_
      - as2.nMoles_/diffNmoles*as2.omega_,
        as1.vibrationalList_,
        as1.nMoles_/diffNmoles*as1.dissociationPotential_
      - as2.nMoles_/diffNmoles*as2.dissociationPotential_,
        as1.noVibrationalTemp_,
        as1.noElectronicLevels_,
        as1.nMoles_/diffNmoles*as1.iHat_
      - as2.nMoles_/diffNmoles*as2.iHat_
    );
}


inline advancedSpecie operator*(const scalar s, const advancedSpecie& as)
{
    return advancedSpecie
    (
        s*as.nMoles_,
        as.molWeight_,
        as.particleType_,
        as.particleCharge_,
        as.diameter_,
        as.omega_,
        as.vibrationalList_,
        as.dissociationPotential_,
        as.noVibrationalTemp_,
        as.noElectronicLevels_,
        as.iHat_
    );
}


inline advancedSpecie operator==(const advancedSpecie& as1, const advancedSpecie& as2)
{
    return as2 - as1;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
