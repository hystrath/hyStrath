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


inline Foam::dsmcParcel::constantProperties::constantProperties
(
    const dictionary& dict
)
:
    mass_(readScalar(dict.lookup("mass"))),
    d_(readScalar(dict.lookup("diameter"))),
    rotationalDegreesOfFreedom_(dict.lookupOrDefault<scalar>("rotationalDegreesOfFreedom", 0)),
    vibrationalDegreesOfFreedom_(dict.lookupOrDefault<scalar>("vibrationalModes", 0)),
    omega_(readScalar(dict.lookup("omega"))),
    alpha_(dict.lookupOrDefault<scalar>("alpha", 1.0)),
    thetaV_(dict.lookupOrDefault<scalarList>("characteristicVibrationalTemperature", scalarList())),
    thetaD_(dict.lookupOrDefault<scalarList>("dissociationTemperature", scalarList())),
    Zref_(dict.lookupOrDefault<scalarList>("Zref", scalarList())),
    TrefZv_(dict.lookupOrDefault<scalarList>("referenceTempForZref", scalarList())),
    charDissQuantumLevel_(dict.lookupOrDefault<labelList>("charDissQuantumLevel", labelList())),
    ionisationTemperature_(dict.lookupOrDefault<scalar>("ionisationTemperature", 1e10)),
    numberOfElectronicLevels_(dict.lookupOrDefault<label>("numberOfElectronicLevels", 1)),
    charge_(dict.lookupOrDefault<label>("charge", 0)),
    electronicEnergyList_(dict.lookupOrDefault<scalarList>("electronicEnergyList", scalarList(1, 0.0))),
    degeneracyList_(dict.lookupOrDefault<labelList>("degeneracyList", labelList(1, 1.0)))
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
    const labelList& vibLevel
)
:
    particle(mesh, position, cellI, tetFaceI, tetPtI),
    U_(U),
    RWF_(RWF),
    ERot_(ERot),
    ELevel_(ELevel),
    typeId_(typeId),
    newParcel_(newParcel),
    //tracked_(nullptr), // TODO TO BE REINTRODUCED
    tracked_(),
    classification_(classification),
    stuck_(nullptr),
    vibLevel_(vibLevel)
{}


// * * * * * * * * * constantProperties Member Functions * * * * * * * * * * //

inline Foam::scalar Foam::dsmcParcel::constantProperties::mass() const
{
    return mass_;
}


inline Foam::scalar
Foam::dsmcParcel::constantProperties::d() const
{
    return d_;
}


inline Foam::scalar Foam::dsmcParcel::constantProperties::sigmaT() const
{
    return constant::mathematical::pi*d_*d_;
}


inline Foam::scalar 
Foam::dsmcParcel::constantProperties::rotationalDegreesOfFreedom() const
{
    return rotationalDegreesOfFreedom_;
}


inline Foam::scalar 
Foam::dsmcParcel::constantProperties::vibrationalDegreesOfFreedom() const
{
    return vibrationalDegreesOfFreedom_;
}


inline Foam::scalar Foam::dsmcParcel::constantProperties::omega() const
{
    return omega_;
}


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

inline const Foam::label& Foam::dsmcParcel::typeId() const
{
    return typeId_;
}


inline const Foam::vector& Foam::dsmcParcel::U() const
{
    return U_;
}


inline const Foam::scalar& Foam::dsmcParcel::RWF() const
{
    return RWF_;
}


inline const Foam::scalar& Foam::dsmcParcel::ERot() const
{
    return ERot_;
}


inline const Foam::label& Foam::dsmcParcel::ELevel() const
{
    return ELevel_;
}


inline const Foam::label& Foam::dsmcParcel::newParcel() const
{
    return newParcel_;
}


inline const Foam::dsmcParcel::TrackedParcel& Foam::dsmcParcel::tracked() const
{
    //return *tracked_; // TODO TO BE REINTRODUCED
    return tracked_;
}


inline bool Foam::dsmcParcel::isTracked() const
{
    /*if(tracked_)
    {
        return true;
    } 
    
    return false;*/// TODO TO BE REINTRODUCED
    
    return tracked_.tracked();
    
}


inline const Foam::label& Foam::dsmcParcel::classification() const
{
    return classification_;
}


inline const Foam::dsmcParcel::StuckParcel& Foam::dsmcParcel::stuck() const
{
    return *stuck_;
}


inline bool Foam::dsmcParcel::isStuck() const
{
    if(stuck_)
    {
        return true;
    }
    
    return false;
}


inline bool Foam::dsmcParcel::isFree() const
{
    if(stuck_)
    {
        return false;
    }
    
    return true;
}


inline const Foam::labelList& Foam::dsmcParcel::vibLevel() const
{
    return vibLevel_;
}


inline Foam::label& Foam::dsmcParcel::typeId()
{
    return typeId_;
}


inline Foam::vector& Foam::dsmcParcel::U()
{
    return U_;
}


inline Foam::scalar& Foam::dsmcParcel::RWF()
{
    return RWF_;
}


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


inline Foam::dsmcParcel::TrackedParcel& Foam::dsmcParcel::tracked()
{
    //return *tracked_; // TODO TO BE REINTRODUCED
    return tracked_;
}


inline void Foam::dsmcParcel::setTracked()
{
    //tracked_ = new dsmcParcel::TrackedParcel(); // TODO TO BE REINTRODUCED
    tracked_ = dsmcParcel::TrackedParcel();
}


inline void Foam::dsmcParcel::setTracked
(
    const bool& tracked,
    const label& inPatchId,
    const scalar& initialTime,
    const vector& initialPosition
)
{
    /*tracked_ = new dsmcParcel::TrackedParcel
        (
            tracked,
            inPatchId,
            initialTime,
            initialPosition
        );*/ // TODO TO BE REINTRODUCED
    
    tracked_ = dsmcParcel::TrackedParcel
        (
            tracked,
            inPatchId,
            initialTime,
            initialPosition
        );    
}


inline void Foam::dsmcParcel::setTracked
(
    const bool& tracked,
    const label& inPatchId,
    const bool& storePositions,
    const scalar& initialTime,
    const vector& initialPosition
)
{
    /*tracked_ = new dsmcParcel::TrackedParcel
        (
            tracked,
            inPatchId,
            storePositions,
            initialTime,
            initialPosition
        );*/ // TODO TO BE REINTRODUCED
        
    tracked_ = dsmcParcel::TrackedParcel
        (
            tracked,
            inPatchId,
            storePositions,
            initialTime,
            initialPosition
        );
}


inline void Foam::dsmcParcel::setTracked
(
    const bool& tracked,
    const label& inPatchId,
    const scalar& initialTime,
    const Foam::vector& initialPosition,
    const Foam::vector& distanceTravelled
)
{
    /*tracked_ = new dsmcParcel::TrackedParcel
        (
            tracked,
            inPatchId,
            initialTime,
            initialPosition,
            distanceTravelled
        );*/ // TODO TO BE REINTRODUCED
        
    tracked_ = dsmcParcel::TrackedParcel
        (
            tracked,
            inPatchId,
            initialTime,
            initialPosition,
            distanceTravelled
        );    
}


inline void Foam::dsmcParcel::setTracked
(
    const bool& tracked,
    const label& inPatchId,
    const bool& storePositions,
    const scalar& initialTime,
    const Foam::vector& initialPosition,
    const Foam::vector& distanceTravelled/*,
    const DynamicList<FixedList<scalar, 4>>& parcelTrajectory*/
)
{
    /*tracked_ = new dsmcParcel::TrackedParcel
        (
            tracked,
            inPatchId,
            storePositions,
            initialTime,
            initialPosition,
            distanceTravelled//,
            //parcelTrajectory
        );*/  // TODO TO BE REINTRODUCED
        
    tracked_ = dsmcParcel::TrackedParcel
        (
            tracked,
            inPatchId,
            storePositions,
            initialTime,
            initialPosition,
            distanceTravelled//,
            //parcelTrajectory
        );    
}


inline void Foam::dsmcParcel::deleteTracked()
{
    //delete tracked_; // TODO TO BE REINTRODUCED
    //tracked_ = nullptr; // TODO TO BE REINTRODUCED
}


inline Foam::label& Foam::dsmcParcel::classification()
{
    return classification_;
}


inline Foam::dsmcParcel::StuckParcel& Foam::dsmcParcel::stuck()
{
    return *stuck_;
}


inline void Foam::dsmcParcel::setStuck()
{
    stuck_ = new dsmcParcel::StuckParcel();
}


inline void Foam::dsmcParcel::setStuck
(
    const Foam::scalarField& wallTemperature,
    const Foam::vectorField& wallVectors
)
{
    stuck_ = new dsmcParcel::StuckParcel
        (
            wallTemperature,
            wallVectors
        );
}


inline void Foam::dsmcParcel::deleteStuck()
{
    delete stuck_;
    stuck_ = nullptr;
}


inline Foam::labelList& Foam::dsmcParcel::vibLevel()
{
    return vibLevel_;
}

// ************************************************************************* //
