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

#include "VariableHardSphere.H"
#include "constants.H"
#include "addToRunTimeSelectionTable.H"

using namespace Foam::constant::mathematical;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(VariableHardSphere, 0);
    addToRunTimeSelectionTable
    (
        BinaryCollisionModel,
        VariableHardSphere,
        dictionary
    );
};


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::VariableHardSphere::VariableHardSphere
(
    const dictionary& dict,
    dsmcCloud& cloud
)
:
    BinaryCollisionModel(dict, cloud),
    coeffDict_
    (
        dict.isDict(typeName + "Coeffs")
        ? dict.subDict(typeName + "Coeffs")
        : dictionary()
    ),
    Tref_(coeffDict_.lookupOrDefault<scalar>("Tref", 273.0))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::VariableHardSphere::~VariableHardSphere()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::VariableHardSphere::active() const
{
    return true;
}


Foam::scalar Foam::VariableHardSphere::sigmaTcR
(
    const dsmcParcel& pP,
    const dsmcParcel& pQ
) const
{
    const label typeIdP = pP.typeId();
    const label typeIdQ = pQ.typeId();

    const scalar dPQ =
        0.5
       *(
            cloud_.constProps(typeIdP).d()
          + cloud_.constProps(typeIdQ).d()
        );

    const scalar omegaPQ =
        0.5
       *(
            cloud_.constProps(typeIdP).omega()
          + cloud_.constProps(typeIdQ).omega()
        );

    const scalar cR = mag(pP.U() - pQ.U());

    if (cR < VSMALL)
    {
        return 0;
    }

    const scalar mP = cloud_.constProps(typeIdP).mass();
    const scalar mQ = cloud_.constProps(typeIdQ).mass();
    const scalar mR = mP*mQ/(mP + mQ);

    //- Calculating cross section = pi*dPQ^2, where dPQ is from Bird, eq. 4.79
    const scalar sigmaTPQ =
        pi*sqr(dPQ)
       *pow(2.0*physicoChemical::k.value()*Tref_/(mR*sqr(cR)), omegaPQ - 0.5)
       /exp(Foam::lgamma(2.5 - omegaPQ));

    return sigmaTPQ*cR;
}


void Foam::VariableHardSphere::collide
(
    dsmcParcel& pP,
    dsmcParcel& pQ,
    const label cellI,
    scalar cR
)
{
    scatter(pP, pQ, cellI, cR);
}


void Foam::VariableHardSphere::scatter
(
    dsmcParcel& pP,
    dsmcParcel& pQ,
    const label cellI,
    scalar cR
)
{
    postCollisionVelocities
    (
        pP.typeId(),
        pQ.typeId(),
        pP.U(),
        pQ.U(),
        cR
    );
    
    //- Collision separation and measurements
    cloud_.cellPropMeasurements().collisionSeparation()[cellI] += 
        sqrt
        (
            sqr(pP.position().x() - pQ.position().x())
          + sqr(pP.position().y() - pQ.position().y())
          + sqr(pP.position().z() - pQ.position().z())
        );
        
    cloud_.cellPropMeasurements().nColls()[cellI]++;
    
    const label classificationP = pP.classification();
    const label classificationQ = pQ.classification();
    
    //- Class I molecule changes to class III molecule when it collides with 
    //  either class II or class III molecules.
    if (classificationP == 0 && classificationQ == 1)
    {
        pP.classification() = 2;
    }
    
    if (classificationQ == 0 && classificationP == 1)
    {
        pQ.classification() = 2;
    }
    
    if (classificationP == 0 && classificationQ == 2)
    {
        pP.classification() = 2;
    }
    
    if (classificationQ == 0 && classificationP == 2)
    {
        pQ.classification() = 2;
    }
}


void Foam::VariableHardSphere::postCollisionVelocities
(
    const label typeIdP,
    const label typeIdQ,
    vector& UP,
    vector& UQ,
    scalar cR
)
{
    if (cR == -1)
    {
        cR = mag(UP - UQ);
    }
    
    const scalar mP = cloud_.constProps(typeIdP).mass();
    const scalar mQ = cloud_.constProps(typeIdQ).mass();

    //- Pre-collision center of mass velocity
    const vector& Ucm = (mP*UP + mQ*UQ)/(mP + mQ);

    const scalar cosTheta = 2.0*cloud_.rndGen().sample01<scalar>() - 1.0;
    const scalar sinTheta = sqrt(1.0 - sqr(cosTheta));
    const scalar phi = twoPi*cloud_.rndGen().sample01<scalar>();

    const vector& postCollisionRelativeU =
        cR
       *vector
        (
            cosTheta,
            sinTheta*cos(phi),
            sinTheta*sin(phi)
        );

    //- Post-collision velocities
    UP = Ucm + postCollisionRelativeU*mQ/(mP + mQ);

    UQ = Ucm - postCollisionRelativeU*mP/(mP + mQ);
}


void Foam::VariableHardSphere::postReactionVelocities
(
    const label typeIdP,
    const label typeIdQ,
    vector& UP,
    vector& UQ,
    scalar cR
)
{
    const scalar mP = cloud_.constProps(typeIdP).mass();
    const scalar mQ = cloud_.constProps(typeIdQ).mass();

    const scalar cosTheta = 2.0*cloud_.rndGen().sample01<scalar>() - 1.0;
    const scalar sinTheta = sqrt(1.0 - sqr(cosTheta));
    const scalar phi = twoPi*cloud_.rndGen().sample01<scalar>();

    const vector& postCollisionRelativeU =
        cR
       *vector
        (
            cosTheta,
            sinTheta*cos(phi),
            sinTheta*sin(phi)
        );

    //- Post-collision velocities
    UQ = UP - postCollisionRelativeU*mP/(mP + mQ);
    
    UP += postCollisionRelativeU*mQ/(mP + mQ);
}


void Foam::VariableHardSphere::redistribute
(
    dsmcParcel& p,
    scalar& translationalEnergy,
    const scalar omegaPQ,
    const bool postReaction 
)
{}


const Foam::dictionary& Foam::VariableHardSphere::coeffDict() const
{
    return coeffDict_;
}

// ************************************************************************* //
