/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2009-2010 OpenCFD Ltd.
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

#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline Foam::dsmcParcel::constantProperties::constantProperties()
:
    mass_(0),
    d_(0)
{}


//template <class ParcelType>
inline Foam::dsmcParcel::constantProperties::constantProperties
(
    const dictionary& dict
)
:
    mass_(readScalar(dict.lookup("mass"))),
    d_(readScalar(dict.lookup("diameter"))),
    rotationalDegreesOfFreedom_(readScalar(dict.lookup("rotationalDegreesOfFreedom"))),
    vibrationalDegreesOfFreedom_(readScalar(dict.lookup("vibrationalModes"))),
    omega_(readScalar(dict.lookup("omega"))),
    alpha_(readScalar(dict.lookup("alpha"))),
    thetaV_(readList<scalar>(dict.lookup("characteristicVibrationalTemperature"))),
    thetaD_(readList<scalar>(dict.lookup("dissociationTemperature"))),
    Zref_(readList<scalar>(dict.lookup("Zref"))),
    TrefZv_(readList<scalar>(dict.lookup("referenceTempForZref"))),
    charDissQuantumLevel_(readList<label>(dict.lookup("charDissQuantumLevel"))),
    ionisationTemperature_(readScalar(dict.lookup("ionisationTemperature"))),
    numberOfElectronicLevels_(readLabel(dict.lookup("numberOfElectronicLevels"))),
    charge_(readLabel(dict.lookup("charge"))),
    electronicEnergyList_(readList<scalar>(dict.lookup("electronicEnergyList"))),
    degeneracyList_(readList<label>(dict.lookup("degeneracyList")))
{
    if(thetaV_.size() != vibrationalDegreesOfFreedom_)
    {
        FatalErrorIn("dsmcParcel::constantProperties::constantProperties")
            << "Number of characteristic vibrational temperatures is "
            << thetaV_.size()
            << ", instead of " 
            << vibrationalDegreesOfFreedom_ << nl 
            << exit(FatalError);
    }
    
    if(thetaD_.size() != vibrationalDegreesOfFreedom_)
    {
        FatalErrorIn("dsmcParcel::constantProperties::constantProperties")
            << "Number of characteristic dissociation temperatures is"
            << thetaD_.size()
            << ", instead of " 
            << vibrationalDegreesOfFreedom_ << nl 
            << exit(FatalError);
    }
    
    if(Zref_.size() != vibrationalDegreesOfFreedom_)
    {
        FatalErrorIn("dsmcParcel::constantProperties::constantProperties")
            << "Number of reference vibrational relaxation numbers is"
            << Zref_.size()
            << ", instead of " 
            << vibrationalDegreesOfFreedom_ << nl 
            << exit(FatalError);
    }
    
    if(TrefZv_.size() != vibrationalDegreesOfFreedom_)
    {
        FatalErrorIn("dsmcParcel::constantProperties::constantProperties")
            << "Number of reference temperature for vibrational relaxation is"
            << TrefZv_.size()
            << ", instead of " 
            << vibrationalDegreesOfFreedom_ << nl 
            << exit(FatalError);
    }
    
    if(charDissQuantumLevel_.size() != vibrationalDegreesOfFreedom_)
    {
        FatalErrorIn("dsmcParcel::constantProperties::constantProperties")
            << "Number of characteristic dissociation quantum levels is"
            << charDissQuantumLevel_.size()
            << ", instead of " 
            << vibrationalDegreesOfFreedom_ << nl 
            << exit(FatalError);
    }
    
    if(degeneracyList_.size() != numberOfElectronicLevels_)
    {
        FatalErrorIn("dsmcParcel::constantProperties::constantProperties")
            << "Number of degeneracy levels should be "
            << numberOfElectronicLevels_
            << ", instead of " 
            << degeneracyList_.size() << nl 
            << exit(FatalError);
    }
    
    if(electronicEnergyList_.size() != numberOfElectronicLevels_)
    {
        FatalErrorIn("dsmcParcel::constantProperties::constantProperties")
            << "Number of electronic energy levels should be "
            << numberOfElectronicLevels_
            << ", instead of " 
            << electronicEnergyList_.size() << nl 
            << exit(FatalError);
    }
    
    if(charge_ > 1 || charge_ < -1)
    {
        FatalErrorIn("dsmcParcel::constantProperties::constantProperties")
            << "Charge value should be 0 for neutrals, "
            << "1 for ions, or -1 for electrons"
            << ", instead of " 
            << charge_ << nl 
            << exit(FatalError);
    }
}



//template <class ParcelType>
inline Foam::dsmcParcel::dsmcParcel
(
    const polyMesh& mesh,
    const vector& position,
    const vector& U,
    const scalar RWF,
    const scalar ERot,
    const label ELevel,
    const label cellI,
    const label tetFaceI,
    const label tetPtI,
    const label typeId,
    const label newParcel,
    const label classification,
    const labelList vibLevel
)
:
    particle(mesh, position, cellI, tetFaceI, tetPtI),
    U_(U),
    RWF_(RWF),
    ERot_(ERot),
    ELevel_(ELevel),
    typeId_(typeId),
    newParcel_(newParcel),
    classification_(classification),
    vibLevel_(vibLevel)
{}


// * * * * * * * * * constantProperties Member Functions * * * * * * * * * * //

//template <class ParcelType>
inline Foam::scalar Foam::dsmcParcel::constantProperties::mass() const
{
    return mass_;
}


//template <class ParcelType>
inline Foam::scalar
Foam::dsmcParcel::constantProperties::d() const
{
    return d_;
}


//template <class ParcelType>
inline Foam::scalar Foam::dsmcParcel::constantProperties::sigmaT() const
{
    return constant::mathematical::pi*d_*d_;
}


//template <class ParcelType>
inline Foam::scalar Foam::dsmcParcel::constantProperties::rotationalDegreesOfFreedom()
const
{
    return rotationalDegreesOfFreedom_;
}

inline Foam::scalar Foam::dsmcParcel::constantProperties::vibrationalDegreesOfFreedom()
const
{
    return vibrationalDegreesOfFreedom_;
}

//template <class ParcelType>
inline Foam::scalar Foam::dsmcParcel::constantProperties::omega() const
{
    return omega_;
}

//template <class ParcelType>
inline Foam::scalar Foam::dsmcParcel::constantProperties::alpha() const
{
    return alpha_;
}

inline Foam::scalarList Foam::dsmcParcel::constantProperties::thetaV() const
{
    return thetaV_;
}

inline Foam::scalarList Foam::dsmcParcel::constantProperties::thetaD() const
{
    return thetaD_;
}

inline Foam::scalarList Foam::dsmcParcel::constantProperties::Zref() const
{
    return Zref_;
}

inline Foam::scalarList Foam::dsmcParcel::constantProperties::TrefZv() const
{
    return TrefZv_;
}

inline Foam::labelList Foam::dsmcParcel::constantProperties::charDissQuantumLevel() const
{
    return charDissQuantumLevel_;
}

inline Foam::scalar Foam::dsmcParcel::constantProperties::ionisationTemperature() const
{
    return ionisationTemperature_;
}

inline Foam::label Foam::dsmcParcel::constantProperties::numberOfElectronicLevels() const
{
    return numberOfElectronicLevels_;
}

inline Foam::label Foam::dsmcParcel::constantProperties::charge() const
{
    return charge_;
}

inline Foam::scalarList Foam::dsmcParcel::constantProperties::electronicEnergyList() const
{
    return electronicEnergyList_;
}

inline Foam::labelList Foam::dsmcParcel::constantProperties::degeneracyList() const
{
    return degeneracyList_;
}


// * * * * * * * * * * * trackData Member Functions  * * * * * * * * * * * * //

//template <class ParcelType>
// inline Foam::dsmcCloud& Foam::dsmcParcel::trackData::cloud()
// {
//     return cloud_;
// }


// * * * * * * * * * * DsmcParcel Member Functions  * * * * * * * * * * //

//template <class ParcelType>
inline Foam::label Foam::dsmcParcel::typeId() const
{
    return typeId_;
}

//template <class ParcelType>
inline const Foam::vector& Foam::dsmcParcel::U() const
{
    return U_;
}

//template <class ParcelType>
inline Foam::scalar Foam::dsmcParcel::RWF() const
{
    return RWF_;
}

//template <class ParcelType>
inline Foam::scalar Foam::dsmcParcel::ERot() const
{
    return ERot_;
}


inline Foam::label Foam::dsmcParcel::ELevel() const
{
    return ELevel_;
}

inline Foam::label Foam::dsmcParcel::newParcel() const
{
    return newParcel_;
}

inline Foam::label Foam::dsmcParcel::classification() const
{
    return classification_;
}

inline Foam::labelList Foam::dsmcParcel::vibLevel() const
{
    return vibLevel_;
}

inline Foam::label& Foam::dsmcParcel::typeId()
{
    return typeId_;
}

//template <class ParcelType>
inline Foam::vector& Foam::dsmcParcel::U()
{
    return U_;
}

//template <class ParcelType>
inline Foam::scalar& Foam::dsmcParcel::RWF()
{
    return RWF_;
}

//template <class ParcelType>
inline Foam::scalar& Foam::dsmcParcel::ERot()
{
    return ERot_;
}

inline Foam::label& Foam::dsmcParcel::ELevel()
{
    return ELevel_;
}

inline Foam::label& Foam::dsmcParcel::newParcel()
{
    return newParcel_;
}

inline Foam::label& Foam::dsmcParcel::classification()
{
    return classification_;
}

inline Foam::labelList& Foam::dsmcParcel::vibLevel()
{
    return vibLevel_;
}
// ************************************************************************* //
