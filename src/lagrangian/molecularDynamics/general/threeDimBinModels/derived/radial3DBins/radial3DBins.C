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
    radial3DBins

Description

\*----------------------------------------------------------------------------*/

#include "radial3DBins.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(radial3DBins, 0);

addToRunTimeSelectionTable(threeDimBinModel, radial3DBins, dictionary);



// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

label radial3DBins::findBinR(const scalar& r)
{
    label n = -1;

    label n2 = label(r/minBinWidth_) - 1;

    if(n2 < 0)
    {
        n2 = 0;
    }

    if(n2 < n_.size())
    {
        n = n_[n2];

        scalar rLimit1 = magRadii_[n] - 0.5*binWidthsR_[n];
        scalar rLimit2 = magRadii_[n] + 0.5*binWidthsR_[n];
    
        if((r >= rLimit1) && (r < rLimit2))
        {}
        else
        {
            n++;
        }
    }

    return n;
}

label radial3DBins::findBinL(const scalar& r)
{
    label n = label(r/binWidthL_);

    if(n == nBinsL_) 
    {
        n--;
    }

    return n;
}

label radial3DBins::findBinA(const scalar& theta)
{
    label n = label(theta/binWidthA_);

    if(n == nBinsA_)
    {
        n--;
    }

    return n;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

//- Construct from components
radial3DBins::radial3DBins
(
    const polyMesh& mesh,
    const dictionary& dict
)
:
    threeDimBinModel(mesh, dict),
    propsDict_(dict.subDict(typeName + "Properties")),
    startPoint_(propsDict_.lookup("startPoint")),
    endPoint_(propsDict_.lookup("endPoint")),
    h_(mag(endPoint_ - startPoint_)),
    unitVector_((endPoint_ - startPoint_)/mag(endPoint_ - startPoint_)),

    angleUnitVectorY_(propsDict_.lookup("angleUnitVectorY")),
    angleUnitVectorX_(propsDict_.lookup("angleUnitVectorX")),
    normalVector_(vector::zero),

    radius_(readScalar(propsDict_.lookup("radius"))),

    nBinsL_(readLabel(propsDict_.lookup("nBinsLength"))),
    nBinsR_(readLabel(propsDict_.lookup("nBinsRadial"))),
    nBinsA_(readLabel(propsDict_.lookup("nBinsAngle"))),

    binWidthL_(h_/scalar(nBinsL_)),
    binWidthR_(radius_/scalar(nBinsR_)),
    binWidthA_(2.0*constant::mathematical::pi/scalar(nBinsA_)),

    magRadii_(),
    binWidthsR_(),
//     volumes_(),
    avVolume_(0.0),
    minBinWidth_(0.0),
    n_()
{

    angleUnitVectorY_ /= mag(angleUnitVectorY_);
    angleUnitVectorX_ /= mag(angleUnitVectorX_);

    //- set volumes --- here we impose a constant volume in all bins,
    // hence the binWidth is going change radially

    avVolume_ = radius_*radius_*constant::mathematical::pi*h_/(nBinsR_*nBinsL_*nBinsA_);

    Info << "average volume: " << avVolume_ << endl;

    DynamicList<scalar> radii(0);
    DynamicList<scalar> binWidthsR(0);

    scalar deltaH = h_/nBinsL_;
    scalar deltaTheta = 360.0/nBinsA_;
    scalar initialBinWidthR = sqrt((avVolume_*360.0)/(deltaTheta*constant::mathematical::pi*deltaH));

    Info << "initialBinWidthR: " << initialBinWidthR << endl;

    binWidthsR.append(initialBinWidthR);

    radii.append(0.5*binWidthsR[0]);

    for(label i = 1; i < nBinsR_; i++)
    {
        scalar r1 = radii[i-1] + 0.5*binWidthsR[i-1];

        scalar r2 = sqrt(((avVolume_*360.0)/(deltaTheta*constant::mathematical::pi*deltaH)) + (r1*r1));

        if(r2 <= radius_)
        {
            radii.append(0.5*(r1+r2));
            binWidthsR.append(r2-r1);
        }
        else
        {
            break;
        }
    }

    //binWidthsR.shrink();
    //radii.shrink();

    nBinsR_ = binWidthsR.size();

    Info << "new nBinsR: " << nBinsR_ << endl;

    Info << "new binWidths: " << binWidthsR << endl;

    magRadii_.setSize(nBinsR_, 0.0);
    binWidthsR_.setSize(nBinsR_, 0.0);
//     volumes_.setSize(nBins_, 0.0);

    forAll(binWidthsR_, n)
    {
        magRadii_[n] = radii[n];
        binWidthsR_[n] = binWidthsR[n];

//         if(n == 0)
//         {
//             volumes_[n] = mathematicalConstant::pi*binWidths_[n]*binWidths_[n]*h_;
//         }
//         else
//         {
//             volumes_[n] = 2.0*mathematicalConstant::pi*magRadii_[n]*binWidths_[n]*h_;
//         }
    }

    Info << "magRadii_: " << magRadii_ << endl;

//     Info << "volumes: " << volumes_ << " avVolume: " << avVolume_ << endl;

    minBinWidth_ = binWidthsR_[nBinsR_-1]/3.0;

    Info << "minBinWidth_: " << minBinWidth_ << " size: " << label(radius_/minBinWidth_) << endl;

    n_.setSize(label(radius_/minBinWidth_), 0);

    forAll(n_, n)
    {
        scalar r = 0.5*minBinWidth_ + scalar(n)*minBinWidth_;

        for(label i = 0; i < nBinsR_; i++)
        {
            scalar rLimit1 = magRadii_[i] - 0.5*binWidthsR_[i];
            scalar rLimit2 = magRadii_[i] + 0.5*binWidthsR_[i];

            if((r >= rLimit1) && (r < rLimit2))
            {
                n_[n] = i;
            }
        }
    }

    normalVector_ = vector(1,0,0);

//     if (propsDict_.found("normalVector"))
//     {
//         normalVector_ = propsDict_.lookup("normalVector");
// 
//         normalVector_ /= mag(normalVector_);
//     }


}



// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

radial3DBins::~radial3DBins()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //






List<label> radial3DBins::isPointWithinBin
(
    const vector& rI,
    const label& cellI
)
{
    List<label> binNumbers;
binNumbers.append(-1);
binNumbers.append(-1);
binNumbers.append(-1);

    vector rSI = rI - startPoint_;
    scalar rD = rSI & unitVector_;

    if((rD <= h_) && (rD >= 0.0))
    {
        scalar rN = mag((rD*unitVector_ + startPoint_) - rI);

        if(rN <= radius_)
        {
            vector rID = rSI - rD*unitVector_;

            scalar theta = acos(angleUnitVectorY_ & rID / mag(rID));

            scalar sign = (angleUnitVectorX_ & rID);

            if(sign < 0.0)
            {
                theta = 2.0*constant::mathematical::pi - theta;
            }

            label nL = findBinL(rD);
            label nR = findBinR(rN);
            label nA = findBinA(theta);

            if
            (
                nR != -1
            )
            {
//                 Pout<< "mol position: " << rI 
//                     << ", rD: " << rD 
//                     << ", bin number " << nL 
//                     << ", radius: " << rN 
//                     << ", bin number: " << nR 
//                     << ", theta (rad): " << theta
//                     << ", theta (deg): " << theta*180.0/mathematicalConstant::pi
//                     << ", bin Number: " << nA
//                     << endl;

                binNumbers[0] = nL;
                binNumbers[1] = nR;
                binNumbers[2] = nA;
            }
        }
    }

    return binNumbers;
}


// scalarField radial3DBins::binPositions()
// {
//     // length
//     scalarField positionsL(nBinsL_, 0.0);
// 
//     forAll(positionsL, i)
//     {
//         positionsL[i] = 0.5*binWidthL_ + scalar(i)*binWidthL_;
//     }
// 
//     // radii
//     scalarField positionsR(nBinsR_, 0.0);
// 
//     forAll(positionsR, i)
//     {
//         positionsR[i] = (0.5 + scalar(i))*binWidthsR_[i];
//     }
// 
//     // angle (radians)
// 
//     scalarField positionsA(nBinsA_, 0.0);
// 
//     forAll(positionsA, i)
//     {
//         positionsA[i] = 0.5*binWidthA_ + scalar(i)*binWidthA_;
//     }
// 
//     return magRadii_;
// }


// 3D positions
vectorField radial3DBins::binPositionsX()
{
    vectorField positionsL(nBinsL_, vector::zero);

    forAll(positionsL, i)
    {
        positionsL[i] = startPoint_ + (0.5 + scalar(i))*binWidthL_*unitVector_;
    }

    return positionsL;
}

scalarField radial3DBins::binPositionsY()
{
    return magRadii_;
}

scalarField radial3DBins::binPositionsZ()
{
    scalarField positionsA(nBinsA_, 0.0);

    forAll(positionsA, i)
    {
        positionsA[i] = (0.5 + scalar(i))*binWidthA_;
    }

    return positionsA;
}


void radial3DBins::write
(
    const fileName& path,
    const word& name
)
{
    vectorField positionsL(nBinsL_, vector::zero);
    const scalarField& positionsR = magRadii_;
    scalarField positionsA(nBinsA_, 0.0);

    forAll(positionsL, i)
    {
        positionsL[i] = startPoint_ + (0.5 + scalar(i))*binWidthL_*unitVector_;
    }

//     forAll(positionsR, i)
//     {
//         positionsR[i] = (0.5 + scalar(i))*binWidthsR_[i];
//     }

    forAll(positionsA, i)
    {
        positionsA[i] = (0.5 + scalar(i))*binWidthA_;
    }

    Info << "path: " << path << endl;

    // output field of positions:

    OFstream positionsFile(path/name+"_positions.xyz");

    label nPositions = nBinsL_ * nBinsR_ * nBinsA_;

    vector h = vector::zero;
    scalar r = 0.0;
    scalar theta = 0.0;
    vector p = vector::zero;

    Info << "radial3DBins::writing out position..." << endl;

    if (positionsFile.good())
    {
        positionsFile << nPositions << endl;
        positionsFile << "(" << endl;

        forAll(positionsL, nL)
        {
            forAll(positionsR, nR)
            {
                forAll(positionsA, nA)
                {
                    h = positionsL[nL];
                    r = positionsR[nR];
                    theta = positionsA[nA];
    
                    p = h + r*cos(theta)*angleUnitVectorY_ 
                        + r*sin(theta)*angleUnitVectorX_;

                    positionsFile 
                        << "(" << p.x() << " " << p.y() << " "
                        << p.z() << ") " << -1
                        << endl;
                }
            }
        }

        positionsFile << ")" << endl;
    }

    Info << "...done." << endl;
}

vector radial3DBins::position(const vector& h, const scalar& r, const scalar& theta)
{
    vector p = h + r*cos(theta)*angleUnitVectorY_ 
                                    + r*sin(theta)*angleUnitVectorX_;
    return p;
}

List<label> radial3DBins::nBins()
{
    List<label> nBins;
    nBins.append(-1);
    nBins.append(-1);
    nBins.append(-1);

    nBins[0] = nBinsL_;
    nBins[1] = nBinsR_;
    nBins[2] = nBinsA_;

    return nBins;
}

scalar radial3DBins::binVolume(const label& n)
{
    return avVolume_;
}


// * * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
