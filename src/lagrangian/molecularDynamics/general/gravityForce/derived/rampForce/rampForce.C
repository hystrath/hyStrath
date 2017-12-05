/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2005 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Class
    rampForce

Description

\*----------------------------------------------------------------------------*/

#include "rampForce.H"
#include "addToRunTimeSelectionTable.H"

// using namespace Foam::constant::mathematical;

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(rampForce, 0);

addToRunTimeSelectionTable(gravityForce, rampForce, dictionary);



// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

//- Construct from components
rampForce::rampForce
(
    Time& time,
    const dictionary& dict
)
:
    gravityForce(time, dict),
    propsDict_(dict.subDict(typeName + "Properties")),
    force_(vector::zero),
    direction_(propsDict_.lookup("forceDirection")),
    initialForce_(readScalar(propsDict_.lookup("initialForce"))),
    finalForce_(readScalar(propsDict_.lookup("finalForce"))),
    initialTime_(readScalar(propsDict_.lookup("initialTime"))),
//     currentTime_(time_.startTime().value()),
    currentTime_(0.0),
    deltaTMD_(time.deltaT().value())
{

    if(propsDict_.found("relaxationTime"))
    {
        relaxationTime_ = readScalar(propsDict_.lookup("relaxationTime"));
    }
    else
    {
        relaxationTime_ = readScalar(propsDict_.lookup("loadingTime"));
    }
//     
    direction_ /= mag(direction_);

    forceGradient_ = (finalForce_ - initialForce_)/relaxationTime_;

    force_ = direction_*initialForce_;
}



// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

rampForce::~rampForce()
{}



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

vector rampForce::force(const vector& position)
{
    return force_;
}

void rampForce::updateForce()
{
//     const scalar t = time_.timeOutputValue();
//     const scalar initialTime = time_.startTime().value();

    currentTime_ += deltaTMD_;

    if(currentTime_ > initialTime_)
    {    
        if(currentTime_ <= relaxationTime_)
        {
            force_ = (forceGradient_*(currentTime_-initialTime_) + initialForce_)*direction_;
        }
        else
        {
                force_ = finalForce_*direction_;
        }
    }

    Info << "force mag = " << mag(force_) << endl;
}

vector rampForce::force(const scalar& time)
{
    return force_;
}

// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //


void rampForce::write
(
    const fileName& fixedPathName,
    const fileName& timePath
)
{

}

void rampForce::updateProperties
(
    const dictionary& dict
)
{}


// * * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
