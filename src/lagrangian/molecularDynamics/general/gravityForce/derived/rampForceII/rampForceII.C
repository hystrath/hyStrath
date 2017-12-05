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
    rampForceII

Description

\*----------------------------------------------------------------------------*/

#include "rampForceII.H"
#include "addToRunTimeSelectionTable.H"

// using namespace Foam::constant::mathematical;

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(rampForceII, 0);

addToRunTimeSelectionTable(gravityForce, rampForceII, dictionary);



// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

//- Construct from components
rampForceII::rampForceII
(
    Time& time,
    const dictionary& dict
)
:
    gravityForce(time, dict),
    propsDict_(dict.subDict(typeName + "Properties")),
    force_(vector::zero),
    direction_(propsDict_.lookup("forceDirection")),
    initialForce_(readScalar(propsDict_.lookup("forceAtTimeZero"))),
    finalForce_(readScalar(propsDict_.lookup("forceAtRampEndTime"))),
    endTime_(readScalar(propsDict_.lookup("rampEndTime"))),
    currentTime_(time.timeOutputValue()),
    deltaTMD_(time.deltaT().value())
{
    timeVarying_ = true;
    
    Info << "current time = " << currentTime_ << endl;

    direction_ /= mag(direction_);

    forceGradient_ = (finalForce_ - initialForce_)/endTime_;
    
    Info << " Force gradient = " << forceGradient_ << endl;
    
    if(propsDict_.found("relaxationTime"))
    {
        relaxationTime_ = readScalar(propsDict_.lookup("relaxationTime"));
        
        bool fixGradient = false;
        
        if (propsDict_.found("fixGradient"))
        {
            fixGradient = Switch(propsDict_.lookup("fixGradient"));
        }
        
        if (relaxationTime_ > 0.0)
        {
            if(!fixGradient)
            {
                forceGradient_ = (finalForce_ - initialForce_)/(endTime_-relaxationTime_);
                
                Info<< nl
                    << " WARNING: Force gradient changed to include initial relaxation time."
                    << " Force gradient = " << forceGradient_ << endl;
            }
            else
            {
                endTime_ += relaxationTime_;
                
                Info<< nl
                    << " WARNING: forceAtRampEndTime changed to include initial relaxation time,"
                    << " = " << endTime_ << endl;
            }
        }
    }
    
    // y intercept
    c_ = finalForce_ - (forceGradient_*endTime_);
    
    force_ = direction_*initialForce_;
    
    
    if(propsDict_.found("acrossMultipleRuns"))
    {
        bool acrossMultipleRuns = Switch(propsDict_.lookup("acrossMultipleRuns"));
        
        if(acrossMultipleRuns)
        {
            force_ = (currentTime_*forceGradient_ + c_)*direction_;

            Info << " ... continuing where we left off ..., force = " << force_ << endl;
        }
    }
    
}



// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

rampForceII::~rampForceII()
{}



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

vector rampForceII::force(const vector& position)
{
    return vector::zero;
}

void rampForceII::updateForce()
{
//     const scalar t = time_.timeOutputValue();
//     const scalar initialTime = time_.startTime().value();

    currentTime_ += deltaTMD_;

    if(currentTime_ > relaxationTime_)
    {    
        if(currentTime_ <= endTime_)
        {
            force_ = ((forceGradient_*currentTime_) + c_)*direction_;
        }
        else
        {
            force_ = finalForce_*direction_;
        }
    }

    Info << "force mag = " << mag(force_) << endl;
}

vector rampForceII::force(const scalar& time)
{
    return force_;
}

// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //


void rampForceII::write
(
    const fileName& fixedPathName,
    const fileName& timePath
)
{

}

void rampForceII::updateProperties
(
    const dictionary& dict
)
{}


// * * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
