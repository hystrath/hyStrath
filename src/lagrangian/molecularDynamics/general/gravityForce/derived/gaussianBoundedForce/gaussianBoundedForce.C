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
    gaussianBoundedForce

Description

\*----------------------------------------------------------------------------*/

#include "gaussianBoundedForce.H"
#include "addToRunTimeSelectionTable.H"
#include "IFstream.H"
#include "writeTimeData.H"


namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(gaussianBoundedForce, 0);

addToRunTimeSelectionTable(gravityForce, gaussianBoundedForce, dictionary);



// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


void gaussianBoundedForce::setBoundBox
(
    const dictionary& propsDict,
    boundedBox& bb,
    const word& name
)
{
    const dictionary& dict(propsDict.subDict(name));
    
    vector startPoint = dict.lookup("startPoint");
    vector endPoint = dict.lookup("endPoint");

    bb.resetBoundedBox(startPoint, endPoint);
    
}




// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

//- Construct from components
gaussianBoundedForce::gaussianBoundedForce
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
//     stress_(readScalar(propsDict_.lookup("stress"))),
    amplitude_(readScalar(propsDict_.lookup("gaussianAmplitude"))),
//     density_(readScalar(propsDict_.lookup("numberDensity"))),
//     constantA_(0.0),
    constantB_(0.0),
    forceDirection_(propsDict_.lookup("forceDirection"))
{
    forceDirection_ /= mag(forceDirection_);

    scalar tF = readScalar(propsDict_.lookup("thickness"));

    scalar sigma = tF/4.0;

//     constantA_ = sqrt(2.0/mathematicalConstant::pi)/sigma;

    constantB_ = 2.0*sigma*sigma;

    setBoundBox(propsDict_, bb_, "forcingRegion");    
    
    // output force distribution

    bool output = false;
    
    if (propsDict_.found("output"))
    {
        output = Switch(propsDict_.lookup("output"));

        if(output)
        {
            const word fieldName = propsDict_.lookup("fieldName");
            scalar nBins = readLabel(propsDict_.lookup("nBins"));
            scalar binWidth = mag(endPoint_ - startPoint_)/(nBins);
        
            scalarField forces(nBins, 0.0);
            scalarField ys(nBins, 0.0);
        
            forAll(forces, n)
            {
                scalar y = binWidth*n;
                ys[n] = y;
//                 forces[n] = ((constantA_*stress_)/density_)*exp(-y*y/constantB_);
                forces[n] = amplitude_*exp(-y*y/constantB_);
            }
        
            // compute integral
            scalar forceIntegral = 0.0;
        
            for (label n = 0; n < forces.size()-1; n++)
            {
                const scalar& r1 = ys[n]; 
                const scalar& r2 = ys[n+1];
        
                const scalar& f1 = forces[n]; 
                const scalar& f2 = forces[n+1];
        
                forceIntegral += 0.5*(r2-r1)*(f1+f2);
            }
        
            Info << "forceIntegral: " << forceIntegral << endl;

            // write out force
            fileName casePath(time.path());
            
            writeTimeData
            (
                casePath,
                "gaussian_"+fieldName+"_force.xy",
                ys,
                forces
            );

            List< Pair<scalar> > forceList(forces.size());

            forAll(forceList, bin)
            {
                forceList[bin].first() = ys[bin];
                forceList[bin].second() = forces[bin];
            }

            writeTimeData
            (
                casePath,
                "gaussian_"+fieldName+"_force",
                forceList
            );
        }
    }
}



// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

gaussianBoundedForce::~gaussianBoundedForce()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

vector gaussianBoundedForce::force(const vector& position)
{
    vector force = vector::zero;

    if(bb_.contains(position))
    {                
        scalar y = (position - startPoint_) & normalVector_;

        force = amplitude_*exp(-y*y/constantB_)*forceDirection_;    
    }

    return force;
}

vector gaussianBoundedForce::force(const scalar& time)
{
    vector force = vector::zero;
    
    return force;
}

void gaussianBoundedForce::updateForce()
{

}

void gaussianBoundedForce::write
(
    const fileName& fixedPathName,
    const fileName& timePath
)
{

}

void gaussianBoundedForce::updateProperties
(
    const dictionary& dict
)
{
    propsDict_ = dict.subDict(typeName + "Properties");

//     stress_ = readScalar(propsDict_.lookup("stress"));
}



// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
