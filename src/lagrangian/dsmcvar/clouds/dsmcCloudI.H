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
inline const Foam::dynamicFvMesh& Foam::dsmcCloud::mesh() const
{
    return mesh_;
}


//template <class dsmcParcel>
inline const Foam::IOdictionary& Foam::dsmcCloud::particleProperties() const
{
    return particleProperties_;
}

//template <class dsmcParcel>
inline const Foam::IOdictionary& Foam::dsmcCloud::controlDict() const
{
    return controlDict_;
}


//template <class dsmcParcel>
inline const Foam::List<Foam::word>& Foam::dsmcCloud::typeIdList() const
{
    return typeIdList_;
}

/*inline Foam::scalar Foam::dsmcCloud::nParticle() const
{
    return nParticle_;
}*/
inline Foam::volScalarField& Foam::dsmcCloud::nParticle()
{
    return nParticle_;
}

inline Foam::scalar Foam::dsmcCloud::nParticle(const label celli)
{
    return nParticle_[celli];
}

inline const Foam::scalarField& Foam::dsmcCloud::deltaTValue()
{
    return variableTimeStep_;
}

inline Foam::scalar Foam::dsmcCloud::deltaTValue(const label celli)
{
    return variableTimeStep_[celli];
}

inline Foam::scalar Foam::dsmcCloud::radialExtent() const
{
    return radialExtent_;
}

inline Foam::scalar Foam::dsmcCloud::maxRWF() const
{
    return maxRWF_;
}

inline bool Foam::dsmcCloud::axisymmetric() const
{
    return axisymmetric_;
}

inline Foam::label Foam::dsmcCloud::nTerminalOutputs() const
{
    return nTerminalOutputs_;
}


//template <class dsmcParcel>
inline const Foam::DynamicList<Foam::DynamicList<Foam::dsmcParcel*> >& Foam::dsmcCloud::cellOccupancy() const // EDITED VINCENT
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

inline const Foam::dsmcControllers& Foam::dsmcCloud::controllers() const  //****
{
    return controllers_;
}

inline Foam::dsmcControllers& Foam::dsmcCloud::controllers()//****
{
    return controllers_;
}

inline const Foam::dsmcDynamicLoadBalancing& Foam::dsmcCloud::dynamicLoadBalancing() const  //****
{
    return dynamicLoadBalancing_;
}

inline Foam::dsmcDynamicLoadBalancing& Foam::dsmcCloud::dynamicLoadBalancing()//****
{
    return dynamicLoadBalancing_;
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
inline Foam::scalarList Foam::dsmcCloud::infoMeasurements() const
{
//     scalar sysMass = 0.0;
    
    scalarList infoMeasurements(5 ,0.0);

    /*forAllConstIter(dsmcCloud, *this, iter)
    {
        const dsmcParcel& p = iter();

        const dsmcParcel::constantProperties& cP = constProps
        (
            p.typeId()
        );

        infoMeasurements[0] += cP.mass(); //massInSystem
        infoMeasurements[1] += 0.5*cP.mass()*(p.U() & p.U()); //linearKineticEnergy
        infoMeasurements[2] += p.ERot(); //rotationalEnergy
        forAll(p.vibLevel(), i)
        {
            infoMeasurements[3] += physicoChemical::k.value()*p.vibLevel()[i]*cP.thetaV()[i]; //vibrationalEnergy
        }
        infoMeasurements[4] += cP.electronicEnergyList()[p.ELevel()]; //electronicEnergy
    }
    
    infoMeasurements[0] *= nParticle_;
    infoMeasurements[1] *= nParticle_;
    infoMeasurements[2] *= nParticle_;
    infoMeasurements[3] *= nParticle_;
    infoMeasurements[4] *= nParticle_;*/
    
    forAllConstIter(dsmcCloud, *this, iter)
    {
        const dsmcParcel& p = iter();

        const dsmcParcel::constantProperties& cP = constProps
        (
            p.typeId()
        );
        
        const scalar nParticle = nParticle_[p.cell()];

        infoMeasurements[0] += cP.mass()*nParticle; //massInSystem
        infoMeasurements[1] += 0.5*cP.mass()*(p.U() & p.U())*nParticle; //linearKineticEnergy
        infoMeasurements[2] += p.ERot()*nParticle; //rotationalEnergy
        forAll(p.vibLevel(), i)
        {
            infoMeasurements[3] += physicoChemical::k.value()*p.vibLevel()[i]*cP.thetaV()[i]*nParticle; //vibrationalEnergy
        }
        infoMeasurements[4] += cP.electronicEnergyList()[p.ELevel()]*nParticle; //electronicEnergy
    }
    
    return infoMeasurements;
}


// //template <class dsmcParcel>
// inline Foam::vector Foam::dsmcCloud::linearMomentumOfSystem() const
// {
//     vector linearMomentum(vector::zero);
// 
//     forAllConstIter(dsmcCloud, *this, iter)
//     {
//         const dsmcParcel& p = iter();
// 
//         const dsmcParcel::constantProperties& cP = constProps
//         (
//             p.typeId()
//         );
// 
//         linearMomentum += cP.mass()*p.U();
//     }
// 
//     return nParticle_*linearMomentum;
// }


// //template <class dsmcParcel>
// inline Foam::scalar Foam::dsmcCloud::linearKineticEnergyOfSystem() const
// {
//     scalar linearKineticEnergy = 0.0;
// 
//     forAllConstIter(dsmcCloud, *this, iter)
//     {
//         const dsmcParcel& p = iter();
// 
//         const dsmcParcel::constantProperties& cP = constProps
//         (
//             p.typeId()
//         );
// 
//         linearKineticEnergy += 0.5*cP.mass()*(p.U() & p.U());
//     }
// 
//     return nParticle_*linearKineticEnergy;
// }
// 
// 
// //template <class dsmcParcel>
// inline Foam::scalar Foam::dsmcCloud::rotationalEnergyOfSystem() const
// {
//     scalar rotationalEnergy = 0.0;
// 
//     forAllConstIter(dsmcCloud, *this, iter)
//     {
//         const dsmcParcel& p = iter();
// 
//         rotationalEnergy += nParticle_*p.ERot();
//     }
// 
//     return rotationalEnergy;
// }
// 
// inline Foam::scalar Foam::dsmcCloud::vibrationalEnergyOfSystem() const
// {
//     scalar vibrationalEnergy = 0.0;
// 
//     forAllConstIter(dsmcCloud, *this, iter)
//     {
//         const dsmcParcel& p = iter();
//         
//         const dsmcParcel::constantProperties& cP = constProps
//         (
//             p.typeId()
//         ); 
// 
//         forAll(p.vibLevel(), i)
//         {
//             vibrationalEnergy += nParticle_*physicoChemical::k.value()*p.vibLevel()[i]*cP.thetaV()[i];
//         }
//     }
// 
//     return vibrationalEnergy;
// }
// 
// inline Foam::scalar Foam::dsmcCloud::electronicEnergyOfSystem() const
// {
//     scalar electronicEnergy = 0.0;
// 
//     forAllConstIter(dsmcCloud, *this, iter)
//     {
//         const dsmcParcel& p = iter();
//         
//         const dsmcParcel::constantProperties& cP = constProps
//         (
//             p.typeId()
//         );
// 
//         electronicEnergy += nParticle_*cP.electronicEnergyList()[p.ELevel()];
//     }
// 
//     return electronicEnergy;
// }

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

inline Foam::cellMeasurements& Foam::dsmcCloud::cellPropMeasurements()
{
    return cellMeas_;
}

inline const Foam::cellMeasurements& Foam::dsmcCloud::cellPropMeasurements() const
{
    return cellMeas_;
}


// ************************************************************************* //
