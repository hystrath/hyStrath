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
    cuboidBinII

Description

\*----------------------------------------------------------------------------*/

#include "cuboidBinII.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(cuboidBinII, 0);

addToRunTimeSelectionTable(binModel, cuboidBinII, dictionary);



// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //
void cuboidBinII::checkBoundBox
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
cuboidBinII::cuboidBinII
(
    const polyMesh& mesh,
    const dictionary& dict
)
:
    binModel(mesh, dict),
    propsDict_(dict.subDict(typeName + "Properties")),
    startPoint_(propsDict_.lookup("startPoint")),
    endPoint_(propsDict_.lookup("endPoint")),
    unitVector_((endPoint_ - startPoint_)/mag(endPoint_ - startPoint_)),
    nBins_(readLabel(propsDict_.lookup("nBins"))),
    binWidth_(mag(endPoint_ - startPoint_)/(nBins_)),

    d1_(readScalar(propsDict_.lookup("d1"))),
    d2_(readScalar(propsDict_.lookup("d2"))),
    n1_(propsDict_.lookup("n1")),
    n2_(propsDict_.lookup("n2")),
    area_(0.0)
{
    n1_ /= mag(n1_);
    n2_ /= mag(n2_);

    area_ = d1_*4.0*d2_;
    
    rSEMag_ = mag(endPoint_ - startPoint_);
    
    midpoint_ = startPoint_ + unitVector_*rSEMag_*0.5;
    
    
    vector minV = startPoint_+ d1_*n1_+d2_*n2_;
    vector maxV = endPoint_ - d1_*n1_ -d2_*n2_;
    
    checkBoundBox(box_, minV, maxV);
}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

cuboidBinII::~cuboidBinII()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


label cuboidBinII::isPointWithinBin
(
    const vector& rI,
    const label& cellI
)
{
    label binNumber = -1;

    if(box_.contains(rI))
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
            
            binNumber = n;            
        }
    }

    return binNumber;
}


scalarField cuboidBinII::binPositions()
{
    scalarField positions(nBins_, 0.0);

    forAll(positions, i)
    {
        positions[i] = 0.5*binWidth_ + scalar(i)*binWidth_;
    }

    return positions;
}

vectorField cuboidBinII::bins()
{
    vectorField positions(nBins_, vector::zero);

    forAll(positions, i)
    {
        positions[i] = startPoint_ + (0.5 + scalar(i))*binWidth_*unitVector_;
    }

    return positions;
}

const label& cuboidBinII::nBins() const
{
    return nBins_;
}

scalar cuboidBinII::binVolume(const label& n)
{
    scalar volume = area_*binWidth_;

    return volume;
}


// * * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
