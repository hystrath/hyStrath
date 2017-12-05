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

#include "constants.H"

using namespace Foam::constant;
using namespace Foam::constant::mathematical;

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


//template <class dsmcParcel>
inline const Foam::word& Foam::dsmcCloud::cloudName() const
{
    return cloudName_;
}


//template <class dsmcParcel>
inline const Foam::fvMesh& Foam::dsmcCloud::mesh() const
{
    return mesh_;
}


//template <class dsmcParcel>
inline const Foam::IOdictionary& Foam::dsmcCloud::particleProperties() const
{
    return particleProperties_;
}


//template <class dsmcParcel>
inline const Foam::List<Foam::word>& Foam::dsmcCloud::typeIdList() const
{
    return typeIdList_;
}


//template <class dsmcParcel>
inline Foam::scalar Foam::dsmcCloud::nParticle() const
{
    return nParticle_;
}


//template <class dsmcParcel>
inline const Foam::List<Foam::DynamicList<Foam::dsmcParcel*> >& Foam::dsmcCloud::cellOccupancy() const
{
    return cellOccupancy_;
}


//template <class dsmcParcel>
inline Foam::volScalarField& Foam::dsmcCloud::sigmaTcRMax()
{
    return sigmaTcRMax_;
}


//template <class dsmcParcel>
inline Foam::scalarField& Foam::dsmcCloud::collisionSelectionRemainder()
{
    return collisionSelectionRemainder_;
}


//template <class dsmcParcel>
inline const Foam::List<Foam::dsmcParcel::constantProperties>& Foam::dsmcCloud::constProps() const
{
    return constProps_;
}


//template <class dsmcParcel>
inline const Foam::dsmcParcel::constantProperties& Foam::dsmcCloud::constProps
(
    label typeId
) const
{
    if (typeId < 0 || typeId >= constProps_.size())
    {
        FatalErrorIn("Foam::dsmcCloud::constProps(label typeId)")
            << "constantProperties for requested typeId index "
            << typeId << " do not exist" << nl
            << abort(FatalError);
    }

    return constProps_[typeId];
}


//template <class dsmcParcel>
inline Foam::Random& Foam::dsmcCloud::rndGen()
{
    return rndGen_;
}

//template <class dsmcParcel>
inline const Foam::BinaryCollisionModel& Foam::dsmcCloud::binaryCollision() const
{
    return binaryCollisionModel_;
}


//template <class dsmcParcel>
inline Foam::BinaryCollisionModel& Foam::dsmcCloud::binaryCollision()
{
    return binaryCollisionModel_();
}

inline const  Foam::dsmcBoundaries&  Foam::dsmcCloud::boundaries() const  //****
{
    return  boundaries_;
}

inline  Foam::dsmcBoundaries&  Foam::dsmcCloud::boundaries()//****
{
    return  boundaries_;
}

inline const Foam::dsmcFieldProperties& Foam::dsmcCloud::fields() const  //****
{
    return fields_;
}

inline Foam::dsmcFieldProperties& Foam::dsmcCloud::fields()//****
{
    return fields_;
}

inline const Foam::dsmcStandardFields& Foam::dsmcCloud::stdFields() const  //****
{
    return standardFields_;
}

inline Foam::dsmcStandardFields& Foam::dsmcCloud::stdFields()//****
{
    return standardFields_;
}

inline const Foam::dsmcControllers& Foam::dsmcCloud::controllers() const  //****
{
    return controllers_;
}

inline Foam::dsmcControllers& Foam::dsmcCloud::controllers()//****
{
    return controllers_;
}

inline const Foam::dsmcFaceTracker& Foam::dsmcCloud::tracker() const//****
{
    return trackingInfo_;
}

inline Foam::dsmcFaceTracker& Foam::dsmcCloud::tracker()//****
{
    return trackingInfo_;
}

inline const Foam::dsmcReactions& Foam::dsmcCloud::reactions() const//****
{
	return reactions_;
}

inline Foam::dsmcReactions& Foam::dsmcCloud::reactions()//****
{
	return reactions_;
}


//template <class dsmcParcel>
inline Foam::scalar Foam::dsmcCloud::massInSystem() const
{
    scalar sysMass = 0.0;

    forAllConstIter(dsmcCloud, *this, iter)
    {
        const dsmcParcel& p = iter();

        const dsmcParcel::constantProperties& cP = constProps
        (
            p.typeId()
        );

        sysMass += cP.mass();
    }

    return nParticle_*sysMass;
}


//template <class dsmcParcel>
inline Foam::vector Foam::dsmcCloud::linearMomentumOfSystem() const
{
    vector linearMomentum(vector::zero);

    forAllConstIter(dsmcCloud, *this, iter)
    {
        const dsmcParcel& p = iter();

        const dsmcParcel::constantProperties& cP = constProps
        (
            p.typeId()
        );

        linearMomentum += cP.mass()*p.U();
    }

    return nParticle_*linearMomentum;
}


//template <class dsmcParcel>
inline Foam::scalar Foam::dsmcCloud::linearKineticEnergyOfSystem() const
{
    scalar linearKineticEnergy = 0.0;

    forAllConstIter(dsmcCloud, *this, iter)
    {
        const dsmcParcel& p = iter();

        const dsmcParcel::constantProperties& cP = constProps
        (
            p.typeId()
        );

        linearKineticEnergy += 0.5*cP.mass()*(p.U() & p.U());
    }

    return nParticle_*linearKineticEnergy;
}


//template <class dsmcParcel>
inline Foam::scalar Foam::dsmcCloud::rotationalEnergyOfSystem() const
{
    scalar rotationalEnergy = 0.0;

    forAllConstIter(dsmcCloud, *this, iter)
    {
        const dsmcParcel& p = iter();

        rotationalEnergy += nParticle_*p.ERot();
    }

    return rotationalEnergy;
}

inline Foam::scalar Foam::dsmcCloud::vibrationalEnergyOfSystem() const
{
    scalar vibrationalEnergy = 0.0;

    forAllConstIter(dsmcCloud, *this, iter)
    {
        const dsmcParcel& p = iter();

        vibrationalEnergy += nParticle_*p.EVib();
    }

    return vibrationalEnergy;
}

//template <class dsmcParcel>
inline Foam::scalar Foam::dsmcCloud::maxwellianAverageSpeed
(
    scalar temperature,
    scalar mass
) const
{
    return
        2.0*sqrt(2.0*physicoChemical::k.value()*temperature/(pi*mass));
}

//template <class dsmcParcel>
inline Foam::scalarField Foam::dsmcCloud::maxwellianAverageSpeed
(
    scalarField temperature,
    scalar mass
) const
{
    return
        2.0*sqrt(2.0*physicoChemical::k.value()*temperature/(pi*mass));
}

//template <class dsmcParcel>
inline Foam::scalar Foam::dsmcCloud::maxwellianRMSSpeed
(
    scalar temperature,
    scalar mass
) const
{
    return sqrt(3.0*physicoChemical::k.value()*temperature/mass);
}

//template <class dsmcParcel>
inline Foam::scalarField Foam::dsmcCloud::maxwellianRMSSpeed
(
    scalarField temperature,
    scalar mass
) const
{
    return sqrt(3.0*physicoChemical::k.value()*temperature/mass);
}


//template <class dsmcParcel>
inline Foam::scalar Foam::dsmcCloud::maxwellianMostProbableSpeed
(
    scalar temperature,
    scalar mass
) const
{
    return sqrt(2.0*physicoChemical::k.value()*temperature/mass);
}


//template <class dsmcParcel>
inline Foam::scalarField Foam::dsmcCloud::maxwellianMostProbableSpeed
(
    scalarField temperature,
    scalar mass
) const
{
    return sqrt(2.0*physicoChemical::k.value()*temperature/mass);
}

//template <class dsmcParcel>
inline void Foam::dsmcCloud::clear()
{
    return IDLList<dsmcParcel>::clear();
}

inline Foam::boundaryMeasurements& Foam::dsmcCloud::boundaryFluxMeasurements()
{
    return boundaryMeas_;
}

inline const Foam::boundaryMeasurements& Foam::dsmcCloud::boundaryFluxMeasurements() const
{
    return boundaryMeas_;
}


// ************************************************************************* //
