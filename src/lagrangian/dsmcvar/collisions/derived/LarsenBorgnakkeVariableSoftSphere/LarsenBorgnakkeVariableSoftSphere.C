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

#include "LarsenBorgnakkeVariableSoftSphere.H"
#include "constants.H"
#include "addToRunTimeSelectionTable.H"

using namespace Foam::constant::mathematical;

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(LarsenBorgnakkeVariableSoftSphere, 0);
    addToRunTimeSelectionTable(BinaryCollisionModel, LarsenBorgnakkeVariableSoftSphere, dictionary);
};

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


Foam::LarsenBorgnakkeVariableSoftSphere::LarsenBorgnakkeVariableSoftSphere
(
    const dictionary& dict,
    dsmcCloud& cloud
)
:
    BinaryCollisionModel(dict, cloud),
    coeffDict_(dict.subDict(typeName + "Coeffs")),
    Tref_(readScalar(coeffDict_.lookup("Tref"))),
    rotationalRelaxationCollisionNumber_
    (
        readScalar(coeffDict_.lookup("rotationalRelaxationCollisionNumber"))
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //


Foam::LarsenBorgnakkeVariableSoftSphere::~LarsenBorgnakkeVariableSoftSphere()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


bool Foam::LarsenBorgnakkeVariableSoftSphere::active() const
{
    return true;
}



Foam::scalar Foam::LarsenBorgnakkeVariableSoftSphere::sigmaTcR
(
    const dsmcParcel& pP,
    const dsmcParcel& pQ
) const
{
//     const CloudType& cloud(this->owner());

    label typeIdP = pP.typeId();
    label typeIdQ = pQ.typeId();

    scalar dPQ =
        0.5
       *(
            cloud_.constProps(typeIdP).d()
          + cloud_.constProps(typeIdQ).d()
        );

    scalar omegaPQ =
        0.5
       *(
            cloud_.constProps(typeIdP).omega()
          + cloud_.constProps(typeIdQ).omega()
        );

    scalar cR = mag(pP.U() - pQ.U());

    if (cR < VSMALL)
    {
        return 0;
    }

    scalar mP = cloud_.constProps(typeIdP).mass();

    scalar mQ = cloud_.constProps(typeIdQ).mass();

    scalar mR = mP*mQ/(mP + mQ);

    // calculating cross section = pi*dPQ^2, where dPQ is from Bird, eq. 4.79
    scalar sigmaTPQ =
        pi*dPQ*dPQ
       *pow(2.0*physicoChemical::k.value()*Tref_/(mR*cR*cR), omegaPQ - 0.5)
       /exp(Foam::lgamma(2.5 - omegaPQ));

    return sigmaTPQ*cR;
}



void Foam::LarsenBorgnakkeVariableSoftSphere::collide
(
    dsmcParcel& pP,
    dsmcParcel& pQ,
    label& cellI
)
{
//     CloudType& cloud_(this->owner());

    label typeIdP = pP.typeId();
    label typeIdQ = pQ.typeId();
    vector& UP = pP.U();
    vector& UQ = pQ.U();
    scalar& ERotP = pP.ERot();
    scalar& ERotQ = pQ.ERot();
    
    scalar alphaPQ = 
    0.5*(
        cloud_.constProps(typeIdP).alpha()
        + cloud_.constProps(typeIdQ).alpha()
    );

    Random& rndGen(cloud_.rndGen());

    scalar inverseCollisionNumber = 1/rotationalRelaxationCollisionNumber_;

    // Larsen Borgnakke rotational energy redistribution part.  Using the serial
    // application of the LB method, as per the INELRS subroutine in Bird's
    // DSMC0R.FOR

    scalar preCollisionERotP = ERotP;

    scalar preCollisionERotQ = ERotQ;

    scalar rotationalDofP = cloud_.constProps(typeIdP).rotationalDegreesOfFreedom();

    scalar rotationalDofQ = cloud_.constProps(typeIdQ).rotationalDegreesOfFreedom();

    scalar omegaPQ =
        0.5
       *(
            cloud_.constProps(typeIdP).omega()
          + cloud_.constProps(typeIdQ).omega()
        );

    scalar mP = cloud_.constProps(typeIdP).mass();

    scalar mQ = cloud_.constProps(typeIdQ).mass();

    scalar mR = mP*mQ/(mP + mQ);

    vector Ucm = (mP*UP + mQ*UQ)/(mP + mQ);

    scalar cRsqr = magSqr(UP - UQ);

    scalar availableEnergy = 0.5*mR*cRsqr;

    scalar ChiB = 2.5 - omegaPQ;

    if (rotationalDofP > 0)
    {
        if (inverseCollisionNumber > rndGen.scalar01())
        {
            availableEnergy += preCollisionERotP;

            scalar ChiA = 0.5*rotationalDofP;

            ERotP = cloud_.energyRatio(ChiA, ChiB)*availableEnergy;

            availableEnergy -= ERotP;
        }
    }

    if (rotationalDofQ > 0)
    {
        if (inverseCollisionNumber > rndGen.scalar01())
        {
            availableEnergy += preCollisionERotQ;
            
            scalar ChiA = 0.5*rotationalDofQ;

            ERotQ = cloud_.energyRatio(ChiA, ChiB)*availableEnergy;

            availableEnergy -= ERotQ;
        }
    }

    // Rescale the translational energy
    scalar A = sqrt(2.0*availableEnergy/mR);
    
    scalar cR = mag(UP - UQ);
    
    vector cRComponents = (UP - UQ) * (A/cR);
    
    cR = A;
    
//     scalar cR = mag(UP - UQ);
    
//     if((alphaPQ - 1.0) < VSMALL)
//     {
//         cR = A;
//     }
//     else
//     {
//         cRComponents *= (A/cR);
//         scalar cR = A;
//     }

    scalar cosTheta = (2.0*pow(rndGen.scalar01(),(1.0/alphaPQ))) - 1.0;

    scalar sinTheta = sqrt(1.0 - cosTheta*cosTheta);

    scalar phi = twoPi*rndGen.scalar01();
    
    scalar D = sqrt(pow(cRComponents.y(),2.0) + pow(cRComponents.z(),2.0));

    vector postCollisionRelU =
        vector
        (
            (cosTheta*cRComponents.x()) + sinTheta*sin(phi)*D,
            (cosTheta*cRComponents.y()) + sinTheta*(cR*cRComponents.z()*cos(phi) - cRComponents.x()*cRComponents.y()*sin(phi))/D,
            (cosTheta*cRComponents.z()) - sinTheta*(cR*cRComponents.y()*cos(phi) + cRComponents.x()*cRComponents.z()*sin(phi))/D
        );

    UP = Ucm + postCollisionRelU*mQ/(mP + mQ);

    UQ = Ucm - postCollisionRelU*mP/(mP + mQ);
}


const Foam::dictionary& Foam::LarsenBorgnakkeVariableSoftSphere::coeffDict() const
{
    return coeffDict_;
}


// ************************************************************************* //
