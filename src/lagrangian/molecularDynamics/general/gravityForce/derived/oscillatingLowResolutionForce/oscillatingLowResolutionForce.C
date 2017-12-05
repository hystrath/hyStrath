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
    oscillatingLowResolutionForce

Description

\*----------------------------------------------------------------------------*/

#include "oscillatingLowResolutionForce.H"
#include "addToRunTimeSelectionTable.H"


// using namespace Foam::constant::mathematical;

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(oscillatingLowResolutionForce, 0);

addToRunTimeSelectionTable(gravityForce, oscillatingLowResolutionForce, dictionary);



// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

//- Construct from components
oscillatingLowResolutionForce::oscillatingLowResolutionForce
(
    Time& time,
    const dictionary& dict
)
:
    gravityForce(time, dict),
    propsDict_(dict.subDict(typeName + "Properties")),
//     force_(propsDict_.lookup("initialForce")),
    force_(vector::zero),
    unitVector_(propsDict_.lookup("unitVector")),
    omega_(readScalar(propsDict_.lookup("omega"))),
    amplitude_(readScalar(propsDict_.lookup("amplitude"))),
    m_(readScalar(propsDict_.lookup("magnitude"))),
//     offsetTime_(0.0),
    currentTime_(time_.startTime().value()),
    deltaTMD_(time.deltaT().value())
{
    unitVector_ /= mag(unitVector_);

//     scalar initialForce = (readScalar(propsDict_.lookup("force")));

//     force_ = unitVector_*initialForce;

//     offsetTime_ = Foam::asin(initialForce/amplitude_)/(360.0*omega_);
}



// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

oscillatingLowResolutionForce::~oscillatingLowResolutionForce()
{}



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

vector oscillatingLowResolutionForce::force(const vector& position)
{
    return force_;
}

void oscillatingLowResolutionForce::updateForce()
{
//     const scalar t = time_.timeOutputValue();
    const scalar initialTime = time_.startTime().value();

    currentTime_ += deltaTMD_;

    scalar time = (currentTime_-initialTime/*+offsetTime_*/);

//     Info << "(floor) label(0.99): " << label(0.99) 
//          << ", (ceil) label(0.99+0.5)" << label(0.99+0.5)
//          << endl;
// 
//     Info << "(floor) label(-0.99): " << label(-0.99)
//          << ", (ceil) label(-0.99-0.5)" << label(-0.99-0.5)
//          << endl;

   // general
    
//     Info << "(floor) label(0.99): " << label(0.99) << ", (ceil) pos: "
//          << label(0.99+(sign(0.99)*0.5)) << ", (ceil) neg: " 
//          << label(-0.99+(sign(-0.99)*0.5))
//          << endl;


//     scalar j = 0.0;

    scalar j = ceil(4.0*m_*((time*omega_)-label(time*omega_)));

//     scalar j= label(v+(sign(v)*0.5));

    force_ = amplitude_*Foam::sin(constant::mathematical::pi*j/(2.0*m_))*unitVector_;

}

vector oscillatingLowResolutionForce::force(const scalar& time)
{
    return force_;
}

void oscillatingLowResolutionForce::write
(
    const fileName& fixedPathName,
    const fileName& timePath
)
{

}


void oscillatingLowResolutionForce::updateProperties
(
    const dictionary& dict
)
{}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
