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
    angularBins

Description

\*----------------------------------------------------------------------------*/

#include "angularBins.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(angularBins, 0);

addToRunTimeSelectionTable(binModel, angularBins, dictionary);



// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

//- Construct from components
angularBins::angularBins
(
    const polyMesh& mesh,
    const dictionary& dict
)
:
    binModel(mesh, dict),
    propsDict_(dict.subDict(typeName + "Properties")),
//     rotationalAxis_(propsDict_.lookup("rotationalAxis")),
    m_(propsDict_.lookup("centrePoint")),
    rVx_(propsDict_.lookup("referenceVectorX")),
    rVy_(propsDict_.lookup("referenceVectorY")),
    Rin_(readScalar(propsDict_.lookup("Rin"))),
    Rout_(readScalar(propsDict_.lookup("Rout"))),
    thetaStart_(readScalar(propsDict_.lookup("angleStart"))),
    thetaEnd_(readScalar(propsDict_.lookup("angleEnd"))),

//     startPoint_(propsDict_.lookup("startPoint")),
//     endPoint_(propsDict_.lookup("endPoint")),
//     unitVector_((endPoint_ - startPoint_)/mag(endPoint_ - startPoint_)),
//     rSEMag_(mag(endPoint_ - startPoint_)),
    nBins_(readLabel(propsDict_.lookup("nBins"))),
//     binWidth_(mag(endPoint_ - startPoint_)/(nBins_)),
//     binWidth_((thetaStart_ - thetaEnd_)/nBins_),
    area_(readScalar(propsDict_.lookup("thickness"))),
    counterClockWise_(false)
{
    scalar deltaR = Rout_ - Rin_;

    area_ *= deltaR;

    if(deltaR <= 0.0)
    {
        FatalErrorIn("angularBins::angularBins()")
            << "Rout " << Rout_ << " needs to be greater than Rin: " << Rin_ 
            << nl << "in: "
            << "system/fieldPropertiesDict"
            << exit(FatalError);
    }

    thetaStart_ *= constant::mathematical::pi/180.0;
    thetaEnd_ *= constant::mathematical::pi/180.0;

    Info<< "thetaStart = " << thetaStart_ 
        << ", thetaEnd = " << thetaEnd_ << endl;

    if(thetaEnd_ <= thetaStart_)
    {
        FatalErrorIn("angularBins::angularBins()")
            << "angleEnd has to be larger than angleStart - define clockwise from referenceVectorY." 
            << nl << "in: "
            << "system/fieldPropertiesDict"
            << exit(FatalError);
    }

    binWidth_ = (thetaEnd_ - thetaStart_)/nBins_;

    Info << "binWidth: " << binWidth_ << endl;

    if (propsDict_.found("counterClockWise"))
    {
        counterClockWise_ = Switch(propsDict_.lookup("counterClockWise"));
    }
}



// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

angularBins::~angularBins()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// cellI is a dummy variable
label angularBins::isPointWithinBin
(
    const vector& rI,
    const label& cellI
)
{
    label binNumber = -1;

    vector rIm = rI - m_;
    vector rImMod = (rVy_ & rIm)*rVy_ + (rVx_ & rIm)*rVx_;

//     Info << "rI = " << rI << ", rIm = " << rIm <<" , rImMod = " << rImMod << endl;

    scalar theta = acos((rVy_ & rImMod)/mag(rImMod));
    scalar sign = rVx_ & rImMod;

    if(sign < 0.0)
    {
        theta = 2.0*constant::mathematical::pi - theta;
    }

//     Info << "theta (deg) = " << theta*180.0/mathematicalConstant::pi << endl;

//     Info << "theta (rad) = " << theta << endl;

    scalar R = mag(rImMod);

//     Info << "R = " << R << endl;

    if ((R >= Rin_) && (R <= Rout_))
    {
        if((theta >= thetaStart_ ) && (theta <= thetaEnd_))
        {
//             Info << "accepted " << endl;

            label n = label((theta-thetaStart_)/binWidth_);

//             Info << "n (test) " << n << endl;    
    
            if(n >= 0) 
            {
                if(n == nBins_) 
                {
                    n--;
                }
        
                if(n < nBins_)
                {
                    binNumber = n;
//                     Info << "bin " << binNumber << endl;
                }
            }
        }
    }

    return binNumber;
}

// angles
scalarField angularBins::binPositions()
{
    scalarField positions(nBins_, 0.0);

    forAll(positions, i)
    {
        positions[i] = 0.5*binWidth_ + scalar(i)*binWidth_;
    }

    return positions;
}

vectorField angularBins::bins()
{
    vectorField positions(nBins_, vector::zero);

    scalar rMean = (Rout_ + Rin_)*0.5;

    forAll(positions, i)
    {
        scalar theta = (0.5 + scalar(i))*binWidth_+thetaStart_;
        positions[i] = m_ + rMean*cos(theta)*rVy_ + rMean*sin(theta)*rVx_;
    }

    return positions;
}

const label& angularBins::nBins() const
{
    return nBins_;
}

scalar angularBins::binVolume(const label& n)
{
    scalar volume = area_*binWidth_*(Rout_ + Rin_)*0.5;

    return volume;
}


// * * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
