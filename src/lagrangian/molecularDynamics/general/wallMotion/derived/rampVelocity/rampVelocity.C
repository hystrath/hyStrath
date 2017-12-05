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
    rampVelocity

Description

\*----------------------------------------------------------------------------*/

#include "rampVelocity.H"
#include "addToRunTimeSelectionTable.H"
#include "mathematicalConstants.H"

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(rampVelocity, 0);

addToRunTimeSelectionTable(wallMotion, rampVelocity, dictionary);



// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

//- Construct from components
rampVelocity::rampVelocity
(
    Time& time,
    const dictionary& dict
)
:
    wallMotion(time, dict),
    propsDict_(dict.subDict(typeName + "Properties")),
    uInitial_(propsDict_.lookup("uInitial")),
    velocity_(uInitial_),
    uMax_(propsDict_.lookup("uMax")),
    tauT_(readScalar(propsDict_.lookup("tauT"))),
    gradient_((uMax_ - uInitial_)/tauT_),
    currentTimeElapsed_(0.0),
    deltaTMD_(time.deltaT().value())
{

}



// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

rampVelocity::~rampVelocity()
{}



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const vector& rampVelocity::velocity() const
{
    return velocity_;
}

void rampVelocity::updateVelocity()
{
//     const scalar t = time_.timeOutputValue();
//     const scalar currentTime = time_.startTime().value();

    currentTimeElapsed_ += deltaTMD_;

//     currentTime_ += deltaTMD_;

//     scalar DeltaT = currentTime_ - initialTime;

    if(currentTimeElapsed_ <= tauT_)
    {
        velocity_ = gradient_*(currentTimeElapsed_) + uInitial_;

        Info << " elapsed time: " << currentTimeElapsed_
             << " velocity: " << velocity_ << endl;

    }
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //
// 
// void rampVelocity::operator=(const rampVelocity& rhs)
// {
//     // Check for assignment to self
//     if (this == &rhs)
//     {
//         FatalErrorIn("rampVelocity::operator=(const rampVelocity&)")
//             << "Attempted assignment to self"
//             << abort(FatalError);
//     }
// 
//     Map<label>::operator=(rhs);
// 
//     binWidth_ = rhs.binWidth();
// }


// * * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

// Ostream& operator<<(Ostream& os, const rampVelocity& d)
// {
//     os  << d.binWidth_
//         << static_cast<const Map<label>&>(d);
// 
//     // Check state of Ostream
//     os.check
//     (
//         "Ostream& operator<<(Ostream&, "
//         "const rampVelocity&)"
//     );
// 
//     return os;
// }

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
