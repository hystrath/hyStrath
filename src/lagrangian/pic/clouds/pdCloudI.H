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

#include "constants.H"

using namespace Foam::constant;
using namespace Foam::constant::mathematical;

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


//template <class pdParcel>
inline const Foam::word& Foam::pdCloud::cloudName() const
{
    return cloudName_;
}


//template <class pdParcel>
inline const Foam::fvMesh& Foam::pdCloud::mesh() const
{
    return mesh_;
}


//template <class pdParcel>
inline const Foam::IOdictionary& Foam::pdCloud::particleProperties() const
{
    return particleProperties_;
}


//template <class pdParcel>
inline const Foam::List<Foam::word>& Foam::pdCloud::typeIdList() const
{
    return typeIdList_;
}


//template <class pdParcel>
inline Foam::scalar Foam::pdCloud::nParticle() const
{
    return nParticle_;
}

inline label Foam::pdCloud::infoFrequency() const
{
    return infoFrequency_;
}

//template <class pdParcel>
inline const Foam::List<Foam::DynamicList<Foam::pdParcel*> >& Foam::pdCloud::cellOccupancy() const
{
    return cellOccupancy_;
}


//template <class pdParcel>
inline Foam::volScalarField& Foam::pdCloud::sigmaTcRMax()
{
    return sigmaTcRMax_;
}

//template <class pdParcel>
inline Foam::scalarField& Foam::pdCloud::collisionSelectionRemainder()
{
    return collisionSelectionRemainder_;
}


//template <class pdParcel>
inline const Foam::List<Foam::pdParcel::constantProperties>& Foam::pdCloud::constProps() const
{
    return constProps_;
}


//template <class pdParcel>
inline const Foam::pdParcel::constantProperties& Foam::pdCloud::constProps
(
    label typeId
) const
{
    if (typeId < 0 || typeId >= constProps_.size())
    {
        FatalErrorIn("Foam::pdCloud::constProps(label typeId)")
            << "constantProperties for requested typeId index "
            << typeId << " do not exist" << nl
            << abort(FatalError);
    }

    return constProps_[typeId];
}


//template <class pdParcel>
inline Foam::Random& Foam::pdCloud::rndGen()
{
    return rndGen_;
}

//template <class pdParcel>
inline const Foam::BinaryCollisionModel& Foam::pdCloud::binaryCollision() const
{
    return binaryCollisionModel_;
}


//template <class pdParcel>
inline Foam::BinaryCollisionModel& Foam::pdCloud::binaryCollision()
{
    return binaryCollisionModel_();
}

inline const  Foam::pdBoundaries&  Foam::pdCloud::boundaries() const  //****
{
    return  boundaries_;
}

inline  Foam::pdBoundaries&  Foam::pdCloud::boundaries()//****
{
    return  boundaries_;
}

inline const Foam::pdFieldProperties& Foam::pdCloud::fields() const  //****
{
    return fields_;
}

inline Foam::pdFieldProperties& Foam::pdCloud::fields()//****
{
    return fields_;
}

inline const Foam::pdStandardFields& Foam::pdCloud::stdFields() const  //****
{
    return standardFields_;
}

inline Foam::pdStandardFields& Foam::pdCloud::stdFields()//****
{
    return standardFields_;
}

inline const Foam::pdEmFields& Foam::pdCloud::emFields() const
{
    return emFields_;
}

inline Foam::pdEmFields& Foam::pdCloud::emFields()
{
    return emFields_;
}

//template <class pdParcel>
inline const Foam::ElectronModel& Foam::pdCloud::electronModel() const
{
    return electronModel_;
}

//template <class pdParcel>
inline Foam::ElectronModel& Foam::pdCloud::electronModel()
{
    return electronModel_();
}
/*
//template <class pdParcel>
inline const Foam::AveragingMethod& Foam::pdCloud::averagingMethod() const
{
    return averagingMethod_;
}

//template <class pdParcel>
inline Foam::AveragingMethod& Foam::pdCloud::averagingMethod()
{
    return averagingMethod_();
}
*/

inline const Foam::pdControllers& Foam::pdCloud::controllers() const  //****
{
    return controllers_;
}

inline Foam::pdControllers& Foam::pdCloud::controllers()//****
{
    return controllers_;
}

inline const Foam::pdFaceTracker& Foam::pdCloud::tracker() const//****
{
    return trackingInfo_;
}

inline Foam::pdFaceTracker& Foam::pdCloud::tracker()//****
{
    return trackingInfo_;
}

inline const Foam::pdReactions& Foam::pdCloud::reactions() const//****
{
	return reactions_;
}

inline Foam::pdReactions& Foam::pdCloud::reactions()//****
{
	return reactions_;
}

/*
//template <class pdParcel>
inline Foam::scalar Foam::pdCloud::massInSystem() const
{
    scalar sysMass = 0.0;

    forAllConstIter(pdCloud, *this, iter)
    {
        const pdParcel& p = iter();

        const pdParcel::constantProperties& cP = constProps
        (
            p.typeId()
        );

        sysMass += cP.mass();
    }

    return nParticle_*sysMass;
}*/
/*
//template <class pdParcel>
inline Foam::scalar Foam::pdCloud::chargeInSystem() const
{
    scalar sysCharge = 0.0;
    const scalar e = electromagnetic::e.value();   //- elementary charge

    forAllConstIter(pdCloud, *this, iter)
    {
        const pdParcel& p = iter();

        const pdParcel::constantProperties& cP = constProps
        (
            p.typeId()
        );

        sysCharge += cP.Ze()*e;
    }

    return nParticle_*sysCharge;
}*/
/*
//template <class pdParcel>
inline Foam::vector Foam::pdCloud::linearMomentumOfSystem() const
{
    vector linearMomentum(vector::zero);

    forAllConstIter(pdCloud, *this, iter)
    {
        const pdParcel& p = iter();

        const pdParcel::constantProperties& cP = constProps
        (
            p.typeId()
        );

        linearMomentum += cP.mass()*p.U();
    }

    return nParticle_*linearMomentum;
}*/

/*
//template <class pdParcel>
inline Foam::scalar Foam::pdCloud::linearKineticEnergyOfSystem() const
{
    scalar linearKineticEnergy = 0.0;

    forAllConstIter(pdCloud, *this, iter)
    {
        const pdParcel& p = iter();

        const pdParcel::constantProperties& cP = constProps
        (
            p.typeId()
        );

        linearKineticEnergy += 0.5*cP.mass()*(p.U() & p.U());
    }

    return nParticle_*linearKineticEnergy;
}*/

/*
//template <class pdParcel>
inline Foam::scalar Foam::pdCloud::rotationalEnergyOfSystem() const
{
    scalar rotationalEnergy = 0.0;

    forAllConstIter(pdCloud, *this, iter)
    {
        const pdParcel& p = iter();

        rotationalEnergy += nParticle_*p.ERot();
    }

    return rotationalEnergy;
}*/
/*
inline Foam::scalar Foam::pdCloud::vibrationalEnergyOfSystem() const
{
    scalar vibrationalEnergy = 0.0;

    forAllConstIter(pdCloud, *this, iter)
    {
        const pdParcel& p = iter();

        vibrationalEnergy += nParticle_*p.EVib();
    }

    return vibrationalEnergy;
}*/

//template <class pdParcel>
inline Foam::scalar Foam::pdCloud::maxwellianAverageSpeed
(
    scalar temperature,
    scalar mass
) const
{
    return
        2.0*sqrt(2.0*physicoChemical::k.value()*temperature/(pi*mass));
}

//template <class pdParcel>
inline Foam::scalarField Foam::pdCloud::maxwellianAverageSpeed
(
    scalarField temperature,
    scalar mass
) const
{
    return
        2.0*sqrt(2.0*physicoChemical::k.value()*temperature/(pi*mass));
}

//template <class pdParcel>
inline Foam::scalar Foam::pdCloud::maxwellianRMSSpeed
(
    scalar temperature,
    scalar mass
) const
{
    return sqrt(3.0*physicoChemical::k.value()*temperature/mass);
}

//template <class pdParcel>
inline Foam::scalarField Foam::pdCloud::maxwellianRMSSpeed
(
    scalarField temperature,
    scalar mass
) const
{
    return sqrt(3.0*physicoChemical::k.value()*temperature/mass);
}


//template <class pdParcel>
inline Foam::scalar Foam::pdCloud::maxwellianMostProbableSpeed
(
    scalar temperature,
    scalar mass
) const
{
    return sqrt(2.0*physicoChemical::k.value()*temperature/mass);
}


//template <class pdParcel>
inline Foam::scalarField Foam::pdCloud::maxwellianMostProbableSpeed
(
    scalarField temperature,
    scalar mass
) const
{
    return sqrt(2.0*physicoChemical::k.value()*temperature/mass);
}

//template <class pdParcel>
inline void Foam::pdCloud::clear()
{
    return IDLList<pdParcel>::clear();
}

inline Foam::boundaryMeasurements& Foam::pdCloud::boundaryFluxMeasurements()
{
    return boundaryMeas_;
}

inline const Foam::boundaryMeasurements& Foam::pdCloud::boundaryFluxMeasurements() const
{
    return boundaryMeas_;
}


// ************************************************************************* //
