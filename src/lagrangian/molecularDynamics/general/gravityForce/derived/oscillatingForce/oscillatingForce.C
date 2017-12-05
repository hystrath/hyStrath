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
    oscillatingForce

Description

\*----------------------------------------------------------------------------*/

#include "oscillatingForce.H"
#include "addToRunTimeSelectionTable.H"

// using namespace Foam::constant::mathematical;

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(oscillatingForce, 0);

addToRunTimeSelectionTable(gravityForce, oscillatingForce, dictionary);



// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

//- Construct from components
oscillatingForce::oscillatingForce
(
    Time& time,
    const dictionary& dict
)
:
    gravityForce(time, dict),
    propsDict_(dict.subDict(typeName + "Properties")),
//     force_(propsDict_.lookup("initialForce")),
    unitVector_(propsDict_.lookup("unitVector")),
    period_(readScalar(propsDict_.lookup("period"))),
    amplitude_(readScalar(propsDict_.lookup("amplitude"))),
    offsetTime_(readScalar(propsDict_.lookup("offsetTime")))
//     currentTime_(time_.startTime().value()),
//     elapsedTime_(0.0)
//     deltaT_(readScalar(propsDict_.lookup("deltaT")))
{
    timeVarying_ = true;
    
    unitVector_ /= mag(unitVector_);

//     scalar initialForce = (readScalar(propsDict_.lookup("force")));

//     force_ = unitVector_*initialForce;

//     scalar t = offsetTime_;
    
//     force_ = amplitude_*Foam::sin(2.0*mathematicalConstant::pi*t/period_)*unitVector_;
}



// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

oscillatingForce::~oscillatingForce()
{}



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

vector oscillatingForce::force(const vector& position)
{
    return vector::zero;
}

void oscillatingForce::updateForce()
{
//     const scalar t = time_.timeOutputValue();
//     const scalar initialTime = time_.startTime().value();

//     elapsedTime_ += deltaT_;
    
//     scalar t = elapsedTime_ + offsetTime_;
    
//     force_ = amplitude_*Foam::sin(2.0*mathematicalConstant::pi*t/period_)*unitVector_;
}



// void oscillatingForce::updateForce(const scalar& time)
// {
// 
// }

vector oscillatingForce::force(const scalar& time)
{
    scalar t = time + offsetTime_;
    
    force_ = amplitude_*Foam::sin(2.0*constant::mathematical::pi*t/period_)*unitVector_;    
    
    return force_;
}

// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //


void oscillatingForce::write
(
    const fileName& fixedPathName,
    const fileName& timePath
)
{

}

void oscillatingForce::updateProperties
(
    const dictionary& dict
)
{}


// * * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
