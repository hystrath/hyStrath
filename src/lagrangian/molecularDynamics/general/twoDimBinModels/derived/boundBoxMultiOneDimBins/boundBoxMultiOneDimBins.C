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
    boundBoxMultiOneDimBins

Description

\*----------------------------------------------------------------------------*/

#include "boundBoxMultiOneDimBins.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(boundBoxMultiOneDimBins, 0);

addToRunTimeSelectionTable(twoDimBinModel, boundBoxMultiOneDimBins, dictionary);



// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void boundBoxMultiOneDimBins::checkBoundBox
(
    boundBox& b,
    const vector& startPoint,
    const vector& endPoint
)
{
    vector& vMin = b.min();
    vector& vMax = b.max();

    if(startPoint.x() < endPoint.x())
    {
        vMin.x() = startPoint.x();
        vMax.x() = endPoint.x();
    }
    else
    {
        vMin.x() = endPoint.x();
        vMax.x() = startPoint.x();
    }
    if(startPoint.y() < endPoint.y())
    {
        vMin.y() = startPoint.y();
        vMax.y() = endPoint.y();
    }
    else
    {
        vMin.y() = endPoint.y();
        vMax.y() = startPoint.y();
    }
    if(startPoint.z() < endPoint.z())
    {
        vMin.z() = startPoint.z();
        vMax.z() = endPoint.z();
    }
    else
    {
        vMin.z() = endPoint.z();
        vMax.z() = startPoint.z();
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

//- Construct from components
boundBoxMultiOneDimBins::boundBoxMultiOneDimBins
(
    const polyMesh& mesh,
    const dictionary& dict
)
:
    twoDimBinModel(mesh, dict),
    propsDict_(dict.subDict(typeName + "Properties")),
    
    // outer bound box
    outerDict_(propsDict_.subDict("outerBoundBox")),
    startPoint_(outerDict_.lookup("startPoint")),
    endPoint_(outerDict_.lookup("endPoint")),

    unitVector_((endPoint_ - startPoint_)/mag(endPoint_ - startPoint_)),
    
    nBins_(readLabel(outerDict_.lookup("nBins"))),
    binWidth_(mag(endPoint_ - startPoint_)/(nBins_)),
    
    d1_(readScalar(outerDict_.lookup("d1"))),
    d2_(readScalar(outerDict_.lookup("d2"))),
    n1_(outerDict_.lookup("n1")),
    n2_(outerDict_.lookup("n2")),
    
    //inner bound boxes
    
    innerDict_(propsDict_.subDict("innerBoundBoxes")),

    nBinsY_(readLabel(innerDict_.lookup("nBins")))
//     unitVectorY_(readLabel(innerDict_.lookup("direction")))

{
    // outer bound box
    
    n1_ /= mag(n1_);
    n2_ /= mag(n2_);

    area_ = d1_*4.0*d2_;    
    
//     rSEMag_ = mag(endPoint_ - startPoint_);
    
    vector minV = startPoint_+ d1_*n1_+d2_*n2_;
    vector maxV = endPoint_ - d1_*n1_ -d2_*n2_;
    
    checkBoundBox(outerBB_, minV, maxV);
    
    vector midPoint = outerBB_.midpoint();
    
    //inner bound boxes
    
    startPointY_ = midPoint - d1_*n1_;
    endPointY_ = midPoint + d1_*n1_;    
    
    Info << "NOTE: direction of measurement of bins is in the n1 direction = " 
          << n1_ << endl;
          
    unitVectorY_ = n1_;
    
    innerBBs_.setSize(nBins_);
    rBins_.setSize(nBins_, vector::zero);
    
    for(label i = 0; i < nBins_; i++)
    {
        rBins_[i] = startPoint_ + unitVector_*binWidth_*0.5 + unitVector_*binWidth_*i;
        
        vector vMin = startPoint_ + unitVector_*binWidth_*i + d1_*n1_+d2_*n2_;
        vector vMax = startPoint_ + unitVector_*binWidth_*(i+1) + d1_*n1_+d2_*n2_;
        
        checkBoundBox(innerBBs_[i], vMin, vMax);        
    }
    
    
    binWidthY_ = (mag(endPointY_ - startPointY_)/(nBinsY_));

    Info << "startpoint Y = " << startPointY_
         << ", endPoint Y = " << endPointY_
         << endl;    
    
    Info << "binWidth X = " << binWidth_
         << ", binWidth Y = " << binWidthY_
         << endl;
}



// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

boundBoxMultiOneDimBins::~boundBoxMultiOneDimBins()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

List<label> boundBoxMultiOneDimBins::isPointWithinBin
(
    const vector& rI,
    const label& cellI
)
{
    List<label> binNumbers;
binNumbers.append(-1);
binNumbers.append(-1);

    bool foundBinX = false;
    bool foundBinY = false;
    
    if(outerBB_.contains(rI))
    {
        {
            vector rSI = rI - startPoint_;    
            scalar rD = rSI & unitVector_;
            label n = label(rD/binWidth_);
            
            if(n >= 0)
            {
                if(n == nBins_) 
                {
                    n--;
                }
                
                binNumbers[0] = n;
                foundBinX = true;
            }
        }
        
        {
            vector rSI = rI - startPointY_;    
            scalar rD = rSI & unitVectorY_;
            label n = label(rD/binWidthY_);
            
            if(n >= 0)
            {
                if(n == nBinsY_) 
                {
                    n--;
                }
                
                binNumbers[1] = n;
                foundBinY = true;
            }        
        }
    }    
    
    if(!foundBinX || !foundBinY)
    {
        binNumbers[0] = -1;
        binNumbers[1] = -1;
    }

    return binNumbers;
}

scalarField boundBoxMultiOneDimBins::binPositionsX()
{
    scalarField positions(nBins_, 0.0);

    forAll(positions, i)
    {
        positions[i] = (0.5 + scalar(i))*binWidth_;
    }
    
    return positions;
}

scalarField boundBoxMultiOneDimBins::binPositionsY()
{
    scalarField positions(nBinsY_, 0.0);

    forAll(positions, i)
    {
        positions[i] = (0.5 + scalar(i))*binWidthY_;
    }

    return positions;
}


void boundBoxMultiOneDimBins::write
(
    const fileName& path,
    const word& name
)
{

}

vector boundBoxMultiOneDimBins::position(/*const vector& h,*/ const scalar& r, const scalar& theta)
{
    vector p = vector::zero;
//     vector p = r*cos(theta)*angleUnitVectorY_ 
//                                     + r*sin(theta)*angleUnitVectorX_;
    return p;
}

const vector& boundBoxMultiOneDimBins::refVector() const
{
    return unitVector_;
}

List<label> boundBoxMultiOneDimBins::nBins()
{
    List<label> nBins;
nBins.append(-1);
nBins.append(-1);

//     nBins[0] = nBinsL_;
    nBins[0] = nBins_;
    nBins[1] = nBinsY_;

    return nBins;
}

scalar boundBoxMultiOneDimBins::binVolume(const label& n)
{
    return binWidth_*binWidthY_*2.0*d2_;
}


// * * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
