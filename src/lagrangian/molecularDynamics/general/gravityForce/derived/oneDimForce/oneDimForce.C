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
    oneDimForce

Description

\*----------------------------------------------------------------------------*/

#include "oneDimForce.H"
#include "addToRunTimeSelectionTable.H"
#include "IFstream.H"


namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(oneDimForce, 0);

addToRunTimeSelectionTable(gravityForce, oneDimForce, dictionary);



// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

//- Construct from components
oneDimForce::oneDimForce
(
    Time& time,
    const dictionary& dict
)
:
    gravityForce(time, dict),
    propsDict_(dict.subDict(typeName + "Properties")),
    startPoint_(propsDict_.lookup("startPoint")),
    endPoint_(propsDict_.lookup("endPoint")),
    normalVector_((endPoint_ - startPoint_)/mag(endPoint_ - startPoint_)),
    forceDirection_(propsDict_.lookup("forceDirection")),
    nBins_(-1),
    binWidth_(-1),
    forces_()

{
    forceDirection_ /= mag(forceDirection_);

    const word distributionName = propsDict_.lookup("forceDistributionName");

//     fileName timePath(time.path()/distributionName);
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
            "void oneDimForce::oneDimForce()"
        )
            << "Cannot open file " << file.name()
            << abort(FatalError);
    }

    nBins_ = forces.size();

    forces_.setSize(nBins_, 0.0);
//     scalar r(nBinsTau_, 0.0);

    forAll(forces, bin)
    {
//         r_[bin] = tau[bin].first();
        forces_[bin] = forces[bin].second();
    }

    binWidth_ = forces[1].first() - forces[0].first();

    length_=mag(startPoint_-endPoint_);

}



// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

oneDimForce::~oneDimForce()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

vector oneDimForce::force(const vector& position)
{
    scalar y = (position - startPoint_) & normalVector_;

    return returnForce(y)*forceDirection_;
}

vector oneDimForce::force(const scalar& time)
{
    vector force = vector::zero;
    
    return force;
}


void oneDimForce::updateForce()
{

}

scalar oneDimForce::returnForce
(
    const scalar& y
)
{
    label n1 = label(y/binWidth_);
    label n2 = n1+1;

    scalar force = 0.0;

    if((n1 >= 0) && (y < length_))
    {
        scalar f1 = forces_[n1];
        scalar f2 = forces_[n2];

        scalar y1 = n1*binWidth_;
        scalar y2 = n2*binWidth_;

        force = f1 + ((f2-f1)*(y-y1)/(y2-y1));
    }

    if(y == length_)
    {
        force = forces_[nBins_-1];
    }

    if(y == 0.0)
    {
        force = forces_[0];
    }
/*
    scalar force = 0.0;
    
    if( (n1 >= 0) && (n2 <= nBins_) )
    {
        if(n1 == nBins_-1)
        {
            force = forces_[n1];
        }
        else if (n2 <= (nBins_-1))
        {
            scalar f1 = forces_[n1];
            scalar f2 = forces_[n2];
            scalar y1 = n1*binWidth_;
            scalar y2 = n2*binWidth_;
    
            force = ((f1 - f2)*(y - y2))/(y1 - y2) + f2;
        }
    }*/

    return force;
}

void oneDimForce::write
(
    const fileName& fixedPathName,
    const fileName& timePath
)
{

}


void oneDimForce::updateProperties
(
    const dictionary& dict
)
{
    propsDict_ = dict.subDict(typeName + "Properties");

    const word distributionName = propsDict_.lookup("forceDistributionName");

//     fileName timePath(time_.path()/distributionName);

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
            "void oneDimForce::oneDimForce()"
        )
            << "Cannot open file " << file.name()
            << abort(FatalError);
    }

    forAll(forces, bin)
    {
//         Info << "original forces: " << forces_[bin] << endl;
        forces_[bin] = forces[bin].second();
//         Info << " new forces: " << forces_[bin] << endl;
    }
}



// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
