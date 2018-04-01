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


inline const Foam::word& Foam::dsmcCloud::cloudName() const
{
    return cloudName_;
}


inline const Foam::dynamicFvMesh& Foam::dsmcCloud::mesh() const
{
    return mesh_;
}


inline const Foam::IOdictionary& Foam::dsmcCloud::particleProperties() const
{
    return particleProperties_;
}


inline const Foam::IOdictionary& Foam::dsmcCloud::controlDict() const
{
    return controlDict_;
}


inline const Foam::List<Foam::word>& Foam::dsmcCloud::typeIdList() const
{
    return typeIdList_;
}


inline Foam::scalar Foam::dsmcCloud::nParticle() const
{
    return nParticle_;
}


inline Foam::scalar Foam::dsmcCloud::radialExtent() const
{
    return radialExtent_;
}


inline Foam::scalar Foam::dsmcCloud::maxRWF() const
{
    return maxRWF_;
}


inline Foam::scalar Foam::dsmcCloud::RWF
(
    const label celli,
    const bool mixedRWM
) const
{
    if (rWMethod_ == "mixed" and mixedRWM)
    {
        return getRWF_cell(celli, mixedRWM);
    }
    
    return RWF_[celli];
}


inline Foam::scalar Foam::dsmcCloud::pRWF
(
    const label patchi,
    const label facei
) const
{
    return RWF_.boundaryField()[patchi][facei];
}


inline bool Foam::dsmcCloud::axisymmetric() const
{
    return axisymmetric_;
}


inline Foam::label Foam::dsmcCloud::nTerminalOutputs() const
{
    return nTerminalOutputs_;
}


inline const Foam::DynamicList<Foam::DynamicList<Foam::dsmcParcel*>>& 
Foam::dsmcCloud::cellOccupancy() const
{
    return cellOccupancy_;
}


inline Foam::volScalarField& Foam::dsmcCloud::sigmaTcRMax()
{
    return sigmaTcRMax_;
}


inline Foam::scalarField& Foam::dsmcCloud::collisionSelectionRemainder()
{
    return collisionSelectionRemainder_;
}


inline const Foam::List<Foam::dsmcParcel::constantProperties>& 
Foam::dsmcCloud::constProps() const
{
    return constProps_;
}


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


inline Foam::Random& Foam::dsmcCloud::rndGen()
{
    return rndGen_;
}


inline const Foam::BinaryCollisionModel& Foam::dsmcCloud::binaryCollision() const
{
    return binaryCollisionModel_;
}


inline Foam::BinaryCollisionModel& Foam::dsmcCloud::binaryCollision()
{
    return binaryCollisionModel_();
}


inline const  Foam::dsmcBoundaries&  Foam::dsmcCloud::boundaries() const
{
    return  boundaries_;
}


inline  Foam::dsmcBoundaries&  Foam::dsmcCloud::boundaries()
{
    return  boundaries_;
}


inline const Foam::dsmcFieldProperties& Foam::dsmcCloud::fields() const
{
    return fields_;
}


inline Foam::dsmcFieldProperties& Foam::dsmcCloud::fields()
{
    return fields_;
}


inline const Foam::dsmcControllers& Foam::dsmcCloud::controllers() const
{
    return controllers_;
}


inline Foam::dsmcControllers& Foam::dsmcCloud::controllers()
{
    return controllers_;
}


inline const Foam::dsmcDynamicLoadBalancing& Foam::dsmcCloud::dynamicLoadBalancing() const
{
    return dynamicLoadBalancing_;
}


inline Foam::dsmcDynamicLoadBalancing& Foam::dsmcCloud::dynamicLoadBalancing()
{
    return dynamicLoadBalancing_;
}


inline const Foam::dsmcFaceTracker& Foam::dsmcCloud::tracker() const
{
    return trackingInfo_;
}


inline Foam::dsmcFaceTracker& Foam::dsmcCloud::tracker()
{
    return trackingInfo_;
}


inline const Foam::dsmcReactions& Foam::dsmcCloud::reactions() const
{
    return reactions_;
}


inline Foam::dsmcReactions& Foam::dsmcCloud::reactions()
{
    return reactions_;
}


inline Foam::scalarList Foam::dsmcCloud::infoMeasurements()
{
    scalarList infoMeasurements(7, 0.0);

    forAllIter(dsmcCloud, *this, iter)
    {
        dsmcParcel& p = iter();
        
        const dsmcParcel::constantProperties& cP = constProps(p.typeId());

        infoMeasurements[0] += cP.mass();
        infoMeasurements[1] += 0.5*cP.mass()*(p.U() & p.U());
        infoMeasurements[2] += p.ERot();
        
        if(cP.vibrationalDegreesOfFreedom() > 0)
        {
            forAll(p.vibLevel(), i)
            {
                infoMeasurements[3] += physicoChemical::k.value()
                    *p.vibLevel()[i]*cP.thetaV()[i];
            }
        }
        
        infoMeasurements[4] += cP.electronicEnergyList()[p.ELevel()];
        
        if(p.isStuck())
        {
            infoMeasurements[5]++;
        }
        else
        {
            infoMeasurements[6]++;
        }
    }
    
    for(int iM=0; iM<=6; iM++)
    {
        infoMeasurements[iM] *= nParticle_;
    }
    
    return infoMeasurements;
}


inline Foam::scalar Foam::dsmcCloud::maxwellianAverageSpeed
(
    const scalar& temperature,
    const scalar& mass
) const
{
    return 2.0*sqrt(2.0*physicoChemical::k.value()*temperature/(pi*mass));
}


inline Foam::scalarField Foam::dsmcCloud::maxwellianAverageSpeed
(
    const scalarField& temperature,
    const scalar& mass
) const
{
    return 2.0*sqrt(2.0*physicoChemical::k.value()*temperature/(pi*mass));
}


inline Foam::scalar Foam::dsmcCloud::maxwellianRMSSpeed
(
    const scalar& temperature,
    const scalar& mass
) const
{
    return sqrt(3.0*physicoChemical::k.value()*temperature/mass);
}


inline Foam::scalarField Foam::dsmcCloud::maxwellianRMSSpeed
(
    const scalarField& temperature,
    const scalar& mass
) const
{
    return sqrt(3.0*physicoChemical::k.value()*temperature/mass);
}


inline Foam::scalar Foam::dsmcCloud::maxwellianMostProbableSpeed
(
    const scalar& temperature,
    const scalar& mass
) const
{
    return sqrt(2.0*physicoChemical::k.value()*temperature/mass);
}


inline Foam::scalarField Foam::dsmcCloud::maxwellianMostProbableSpeed
(
    const scalarField& temperature,
    const scalar& mass
) const
{
    return sqrt(2.0*physicoChemical::k.value()*temperature/mass);
}


inline void Foam::dsmcCloud::clear()
{
    return IDLList<dsmcParcel>::clear();
}


inline Foam::boundaryMeasurements& Foam::dsmcCloud::boundaryFluxMeasurements()
{
    return boundaryMeas_;
}


inline const Foam::boundaryMeasurements& 
Foam::dsmcCloud::boundaryFluxMeasurements() const
{
    return boundaryMeas_;
}


inline Foam::cellMeasurements& Foam::dsmcCloud::cellPropMeasurements()
{
    return cellMeas_;
}


inline const Foam::cellMeasurements& 
Foam::dsmcCloud::cellPropMeasurements() const
{
    return cellMeas_;
}


inline bool Foam::dsmcCloud::measureMediumProperties() const
{
    return measureMediumProperties_;
}


inline bool Foam::dsmcCloud::measureInstantaneousMSD() const
{
    return measureInstantaneousMSD_;
}


// ************************************************************************* //
