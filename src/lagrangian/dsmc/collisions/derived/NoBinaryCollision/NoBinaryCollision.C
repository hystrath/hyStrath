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

#include "NoBinaryCollision.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(NoBinaryCollision, 0);
    addToRunTimeSelectionTable
    (
        BinaryCollisionModel,
        NoBinaryCollision,
        dictionary
    );
};

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::NoBinaryCollision::NoBinaryCollision
(
    const dictionary& dict,
    dsmcCloud& cloud
)
:
    BinaryCollisionModel(cloud)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::NoBinaryCollision::~NoBinaryCollision()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


bool Foam::NoBinaryCollision::active() const
{
    return false;
}


Foam::scalar Foam::NoBinaryCollision::sigmaTcR
(
    const dsmcParcel& pP,
    const dsmcParcel& pQ
) const
{
    FatalErrorIn
    (
        "Foam::scalar Foam::NoBinaryCollision::sigmaTcR"
        "("
            "const dsmcParcel& pP,"
            "const dsmcParcel& pQ"
        ") const"
    )
        << "sigmaTcR called in NoBinaryCollision model, this should "
        << "not happen, this is not an actual model." << nl
        << "Enclose calls to sigmaTcR within a if (binaryCollision().active()) "
        << "check."
        << abort(FatalError);

    return 0.0;
}


void Foam::NoBinaryCollision::collide
(
    dsmcParcel& pP,
    dsmcParcel& pQ,
    const label cellI,
    scalar cR
)
{}


void Foam::NoBinaryCollision::scatter
(
    dsmcParcel& pP,
    dsmcParcel& pQ,
    const label cellI,
    scalar cR
)
{}


void Foam::NoBinaryCollision::redistribute
(
    dsmcParcel& p,
    scalar& translationalEnergy,
    const scalar omegaPQ,
    const bool postReaction
)
{}


const Foam::dictionary& Foam::NoBinaryCollision::coeffDict() const
{
    return coeffDict_;
}


// ************************************************************************* //
