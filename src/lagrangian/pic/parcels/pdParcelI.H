/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016-2021 hyStrath
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

#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline Foam::pdParcel::constantProperties::constantProperties()
:
    mass_(0),
    Ze_(0),
    d_(0)
{}

/*inline Foam::pdParcel::constantProperties::~constantProperties()
{
    //Info << "pdParcel::constantProperties Destructor" << endl;
}*/

//template <class ParcelType>
inline Foam::pdParcel::constantProperties::constantProperties
(
    const dictionary& dict
)
:
    mass_(readScalar(dict.lookup("mass"))),
    Ze_(readScalar(dict.lookup("Ze"))),
    d_(readScalar(dict.lookup("diameter"))),
    rotationalDegreesOfFreedom_
    (
        dict.lookupOrDefault<scalar>("rotationalDegreesOfFreedom", 0)
    ),
    vibrationalDegreesOfFreedom_
    (
        dict.lookupOrDefault<scalar>("vibrationalModes", 0)
    ),
    omega_(readScalar(dict.lookup("omega"))),
    alpha_(dict.lookupOrDefault<scalar>("alpha", 1.0)),
    thetaV_(dict.lookupOrDefault<scalar>("characteristicVibrationalTemperature", 0.0)),
    charDissQuantumLevel_(dict.lookupOrDefault<label>("charDissQuantumLevel", 0.0)),
    thetaD_(dict.lookupOrDefault<scalar>("dissociationTemperature", 0.0)),
    Zref_(dict.lookupOrDefault<scalar>("Zref", 0.0)),
    TrefZv_(dict.lookupOrDefault<scalar>("referenceTempForZref", 0.0))
{}

inline Foam::pdParcel::trackingData::trackingData(pdCloud& cloud, label part)
:
    particle::TrackingData<pdCloud>(cloud),
    part_(part)
{}

//template <class ParcelType>
inline Foam::pdParcel::pdParcel
(
    const polyMesh& mesh,
    const vector& position,
    const vector& U,
    const vector& A,
    const scalar EPot,
    const scalar ERot,
    const scalar EVib,
    const label cellI,
    const label tetFaceI,
    const label tetPtI,
    const label typeId,
    const label newParcel,
    const label classification
)
:
    particle(mesh, position, cellI, tetFaceI, tetPtI),
    U_(U),
    A_(A),
    EPot_(EPot),
    ERot_(ERot),
    EVib_(EVib),
    typeId_(typeId),
    newParcel_(newParcel),
    classification_(classification)
{}


// * * * * * * * * * constantProperties Member Functions * * * * * * * * * * //

//template <class ParcelType>
inline Foam::scalar Foam::pdParcel::constantProperties::mass() const
{
    return mass_;
}

//template <class ParcelType>
inline Foam::scalar Foam::pdParcel::constantProperties::Ze() const
{
    return Ze_;
}

//template <class ParcelType>
inline Foam::scalar
Foam::pdParcel::constantProperties::d() const
{
    return d_;
}

//template <class ParcelType>
inline Foam::scalar Foam::pdParcel::constantProperties::sigmaT() const
{
    return constant::mathematical::pi*d_*d_;
}

//template <class ParcelType>
inline Foam::scalar Foam::pdParcel::constantProperties::rotationalDegreesOfFreedom()
const
{
    return rotationalDegreesOfFreedom_;
}

inline Foam::scalar Foam::pdParcel::constantProperties::vibrationalDegreesOfFreedom()
const
{
    return vibrationalDegreesOfFreedom_;
}


inline Foam::label Foam::pdParcel::constantProperties::charDissQuantumLevel() const
{
    return charDissQuantumLevel_;
}

//template <class ParcelType>
inline Foam::scalar Foam::pdParcel::constantProperties::omega() const
{
    return omega_;
}

//template <class ParcelType>
inline Foam::scalar Foam::pdParcel::constantProperties::alpha() const
{
    return alpha_;
}

inline Foam::scalar Foam::pdParcel::constantProperties::thetaV() const
{
    return thetaV_;
}

inline Foam::scalar Foam::pdParcel::constantProperties::thetaD() const
{
    return thetaD_;
}

inline Foam::scalar Foam::pdParcel::constantProperties::Zref() const
{
    return Zref_;
}

inline Foam::scalar Foam::pdParcel::constantProperties::TrefZv() const
{
    return TrefZv_;
}


// * * * * * * * * * * * trackData Member Functions  * * * * * * * * * * * * //

inline Foam::label Foam::pdParcel::trackingData::part() const
{
    return part_;
}

// * * * * * * * * * * PdParcel Member Functions  * * * * * * * * * * //

//template <class ParcelType>
inline Foam::label Foam::pdParcel::typeId() const
{
    return typeId_;
}

//template <class ParcelType>
inline const Foam::vector& Foam::pdParcel::U() const
{
    return U_;
}

//template <class ParcelType>
inline const Foam::vector& Foam::pdParcel::A() const
{
    return A_;
}

//template <class ParcelType>
inline Foam::scalar Foam::pdParcel::EPot() const
{
    return EPot_;
}

//template <class ParcelType>
inline Foam::scalar Foam::pdParcel::ERot() const
{
    return ERot_;
}

inline Foam::scalar Foam::pdParcel::EVib() const
{
    return EVib_;
}

inline Foam::label Foam::pdParcel::newParcel() const
{
    return newParcel_;
}

inline Foam::label Foam::pdParcel::classification() const
{
    return classification_;
}

inline Foam::label& Foam::pdParcel::typeId()
{
    return typeId_;
}

//template <class ParcelType>
inline Foam::vector& Foam::pdParcel::U()
{
    return U_;
}

inline Foam::vector& Foam::pdParcel::A()
{
    return A_;
}

//template <class ParcelType>
inline Foam::scalar& Foam::pdParcel::EPot()
{
    return EPot_;
}

//template <class ParcelType>
inline Foam::scalar& Foam::pdParcel::ERot()
{
    return ERot_;
}

inline Foam::scalar& Foam::pdParcel::EVib()
{
    return EVib_;
}

inline Foam::label& Foam::pdParcel::newParcel()
{
    return newParcel_;
}

inline Foam::label& Foam::pdParcel::classification()
{
    return classification_;
}
