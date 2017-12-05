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
    timeVaryingForce

Description

\*----------------------------------------------------------------------------*/

#include "timeVaryingForce.H"
#include "addToRunTimeSelectionTable.H"
#include "IFstream.H"


namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(timeVaryingForce, 0);

addToRunTimeSelectionTable(gravityForce, timeVaryingForce, dictionary);



// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

//- Construct from components
timeVaryingForce::timeVaryingForce
(
    Time& time,
    const dictionary& dict
)
:
    gravityForce(time, dict),
    propsDict_(dict.subDict(typeName + "Properties")),
    force_(vector::zero),
    forceDirection_(propsDict_.lookup("forceDirection")),
    offsetTime_(readScalar(propsDict_.lookup("offsetTime"))),
//     deltaT_(readScalar(propsDict_.lookup("deltaT"))),    
//     elapsedTime_(0.0),
//     index_(0),
//     times_(),
    forces_()

        
{
    timeVarying_ = true;
    
    forceDirection_ /= mag(forceDirection_);

    const word distributionName = propsDict_.lookup("forceDistributionName");

    fileName timePath(time_.time().system()/distributionName);

    IFstream file(timePath);

    List< Pair<scalar> > forces;

    if (file.good())
    {
        file >> forces;
    }
    else
    {
        FatalErrorIn
        (
            "void timeVaryingForce::timeVaryingForce()"
        )
            << "Cannot open file " << file.name()
            << abort(FatalError);
    }

    nBins_ = forces.size();
    
    forces_.setSize(nBins_, 0.0);

    forAll(forces, bin)
    {
        forces_[bin] = forces[bin].second();
    }
    
    binWidth_ = forces[1].first()-forces[0].first();

    Info << "binWidth = " << binWidth_ << endl;
}



// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

timeVaryingForce::~timeVaryingForce()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

vector timeVaryingForce::force(const vector& position)
{
    return vector::zero;
}

vector timeVaryingForce::force(const scalar& time)
{
    scalar t = time + offsetTime_;
    
    force_ = getForce(t);
    
//     Info << "force = " << force_ << endl;
    
    return force_;
}

void timeVaryingForce::updateForce()
{}

vector timeVaryingForce::getForce(const scalar& t)
{
    label index = label(t/binWidth_);
    
    if(index < nBins_)
    {
        return forces_[index]*forceDirection_;
    }
    else
    {
        Info<< "WARNING in timeVaryingForce::getForce() " 
            <<  nl << "exceeded the time-varying list. "
            <<endl;
            
        return forces_[nBins_-1]*forceDirection_;    
    }
}

void timeVaryingForce::write
(
    const fileName& fixedPathName,
    const fileName& timePath
)
{

}


void timeVaryingForce::updateProperties
(
    const dictionary& dict
)
{
    propsDict_ = dict.subDict(typeName + "Properties");
}



// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
