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
    oscillatingVelocity

Description

\*----------------------------------------------------------------------------*/

#include "oscillatingVelocity.H"
#include "addToRunTimeSelectionTable.H"
#include "mathematicalConstants.H"

// using namespace Foam::constant::mathematical;

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(oscillatingVelocity, 0);

addToRunTimeSelectionTable(wallMotion, oscillatingVelocity, dictionary);



// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

//- Construct from components
oscillatingVelocity::oscillatingVelocity
(
    Time& time,
    const dictionary& dict
)
:
    wallMotion(time, dict),
    propsDict_(dict.subDict(typeName + "Properties")),
    velocity_(vector::zero),
    unitVector_(propsDict_.lookup("unitVector")),    
    initialVelocity_(readScalar(propsDict_.lookup("initialVelocity"))),
    amplitude_(readScalar(propsDict_.lookup("amplitude"))),
    period_(readScalar(propsDict_.lookup("period"))),
    offsetTime_(0.0),
    currentTime_(time_.startTime().value()),
//     deltaTMD_(time.deltaT().value()),
    deltaT_(readScalar(propsDict_.lookup("deltaT")))    
{
    unitVector_ /= mag(unitVector_);
    
    velocity_ = unitVector_*initialVelocity_;    
    
//     offsetTime_ = tauT_*Foam::asin(mag(velocity_)/mag(uMax_))/360;
}



// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

oscillatingVelocity::~oscillatingVelocity()
{}



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const vector& oscillatingVelocity::velocity() const
{
    return velocity_;
}

void oscillatingVelocity::updateVelocity()
{
//     const scalar t = time_.timeOutputValue();
    const scalar initialTime = time_.startTime().value();
    
    currentTime_ += deltaT_;
    
    scalar t = currentTime_-initialTime+offsetTime_;
    
    velocity_ = amplitude_*Foam::sin(2.0*constant::mathematical::pi*t/period_)*unitVector_;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

// Ostream& operator<<(Ostream& os, const oscillatingVelocity& d)
// {
//     os  << d.binWidth_
//         << static_cast<const Map<label>&>(d);
// 
//     // Check state of Ostream
//     os.check
//     (
//         "Ostream& operator<<(Ostream&, "
//         "const oscillatingVelocity&)"
//     );
// 
//     return os;
// }

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
