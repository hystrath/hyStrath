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
    mixedOscillatingForce

Description

\*----------------------------------------------------------------------------*/

#include "mixedOscillatingForce.H"
#include "addToRunTimeSelectionTable.H"

// using namespace Foam::constant::mathematical;

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(mixedOscillatingForce, 0);

addToRunTimeSelectionTable(gravityForce, mixedOscillatingForce, dictionary);



// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

//- Construct from components
mixedOscillatingForce::mixedOscillatingForce
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
//     period_(readScalar(propsDict_.lookup("period"))),
    amplitude_(readScalar(propsDict_.lookup("amplitude"))),
//     offsetTime_(0.0),
    currentTime_(0.0),
//     currentTime_(time_.startTime().value()),
    deltaT_(readScalar(propsDict_.lookup("deltaT")))
{
    unitVector_ /= mag(unitVector_);
        
    periodCoeffs_ = List<scalar>(propsDict_.lookup("coeffs"));   
    
//     scalar initialForce = (readScalar(propsDict_.lookup("force")));

//     force_ = unitVector_*initialForce;

    //offsetTime_ = Foam::asin(initialForce/amplitude_)/(360.0*omega_);
    
    
    bool outputForces = false;
    
    if (propsDict_.found("output"))
    {
        outputForces = Switch(propsDict_.lookup("output"));    
        
        if(outputForces)
        {    
            output(time);
        }
    }
}



// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

mixedOscillatingForce::~mixedOscillatingForce()
{}



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

vector mixedOscillatingForce::force(const vector& position)
{
    return force_;
}

void mixedOscillatingForce::updateForce()
{
//     const scalar t = time_.timeOutputValue();
//     const scalar initialTime = time_.startTime().value();

    currentTime_ += deltaT_;
    
//     scalar t = currentTime_-initialTime+offsetTime_;
    
    scalar t = currentTime_;
    
    scalar period = getPeriod(t);
    
    force_ = amplitude_*Foam::sin(2.0*constant::mathematical::pi*t/period)*unitVector_;
}

vector mixedOscillatingForce::force(const scalar& time)
{
    return force_;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //
scalar mixedOscillatingForce::getPeriod(const scalar& t)
{
    scalar period = 0.0;
    
    label M = periodCoeffs_.size() - 1;
    
    forAll(periodCoeffs_, j)
    {
        period += periodCoeffs_[j]*pow(t, (M-j));
    }
    
    return period;
}

void mixedOscillatingForce::output(Time& time)
{
    label N = 4e6;
    
    scalarField t(N, 0.0);     
    scalarField period(N, 0.0);   
    scalarField forces(N, 0.0);
    
    for (label i=0; i< N; i++)
    {
        updateForce();
        
        t[i] = currentTime_;
        forces[i]=force_.x();
        
        scalar p = getPeriod(currentTime_);
        period[i]=p;        
    }
    
    fileName casePath(time.path());   

    writeTimeData
    (
        casePath,
        "mixedOscillatingForce_forces.xy",
        t,
        forces
    );

    writeTimeData
    (
        casePath,
        "mixedOscillatingForce_period.xy",
        t,
        period
    );    
    
    
}

void mixedOscillatingForce::write
(
    const fileName& fixedPathName,
    const fileName& timePath
)
{

}

void mixedOscillatingForce::updateProperties
(
    const dictionary& dict
)
{}


// * * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
