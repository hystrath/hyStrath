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
    multiRampForce

Description

\*----------------------------------------------------------------------------*/

#include "multiRampForce.H"
#include "addToRunTimeSelectionTable.H"

// using namespace Foam::constant::mathematical;

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(multiRampForce, 0);

addToRunTimeSelectionTable(gravityForce, multiRampForce, dictionary);



// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

//- Construct from components
multiRampForce::multiRampForce
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
    loadingTime_(readScalar(propsDict_.lookup("loadingTime"))),
    breathingTime_(readScalar(propsDict_.lookup("breathingTime"))),
    noOfStops_(readLabel(propsDict_.lookup("noOfStops"))),
    currentTime_(0.0),
    timeLoading_(0.0),
    timeBreathing_(0.0),
    deltaTMD_(time.deltaT().value())

{
    partialLoadingTime_ = loadingTime_/noOfStops_;
    
    direction_ /= mag(direction_);

    forceGradient_ = (finalForce_ - initialForce_)/loadingTime_;

    force_ = direction_*initialForce_;
    
    Info << nl << "breathing time per stop = " << breathingTime_
         << ", total breathing time = " << breathingTime_*noOfStops_
         << nl << "loading time per stop = " << partialLoadingTime_
         << ", total loading time = " << loadingTime_
         << nl << "Total time = " << (breathingTime_*noOfStops_) + loadingTime_
         << endl;
}



// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

multiRampForce::~multiRampForce()
{}



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

vector multiRampForce::force(const vector& position)
{
    return force_;
}

void multiRampForce::updateForce()
{
//     const scalar t = time_.timeOutputValue();
//     const scalar initialTime = time_.startTime().value();

    currentTime_ += deltaTMD_;
    timeLoading_ += deltaTMD_;
    
    if(currentTime_ <= loadingTime_)
    {
        if(timeLoading_ <= partialLoadingTime_)
        {
            Info << "loading" << endl;
            
            force_ = (forceGradient_*currentTime_ + initialForce_)*direction_;
            timeBreathing_ = 0.0;
        }
        else
        {
            Info << "breathing" << endl;
            
            timeBreathing_ += deltaTMD_;
            
            initialForce_ -= forceGradient_*deltaTMD_;
            
            if(timeBreathing_ >= breathingTime_)
            {
                timeLoading_ = 0.0;
            }
        }
    }
    else
    {
        force_ = finalForce_*direction_;
    }
    
    Info << "force mag = " << mag(force_) << endl;
}



vector multiRampForce::force(const scalar& time)
{
    return force_;
}
// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //


void multiRampForce::write
(
    const fileName& fixedPathName,
    const fileName& timePath
)
{

}

void multiRampForce::updateProperties
(
    const dictionary& dict
)
{}


// * * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
