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
    twoDimBins

Description

\*----------------------------------------------------------------------------*/

#include "twoDimBins.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(twoDimBins, 0);

addToRunTimeSelectionTable(twoDimBinModel, twoDimBins, dictionary);



// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

//- Construct from components
twoDimBins::twoDimBins
(
    const polyMesh& mesh,
    const dictionary& dict
)
:
    twoDimBinModel(mesh, dict),
    propsDict_(dict.subDict(typeName + "Properties")),
    startPoint_(propsDict_.lookup("startPoint")),
    endPoint_(propsDict_.lookup("endPoint")),
    h_(mag(endPoint_ - startPoint_)),
    unitVector_((endPoint_ - startPoint_)/mag(endPoint_ - startPoint_)),

    unitVectorX_(propsDict_.lookup("unitVectorX")),
    unitVectorY_(propsDict_.lookup("unitVectorY")),

    nBinsX_(readLabel(propsDict_.lookup("nBinsX"))),
    nBinsY_(readLabel(propsDict_.lookup("nBinsY"))),

    dX_(readScalar(propsDict_.lookup("dX"))),
    dY_(readScalar(propsDict_.lookup("dY"))),
    
    binWidthX_(2.0*dX_/nBinsX_),
    binWidthY_(2.0*dY_/nBinsY_)

{

    unitVectorX_ /= mag(unitVectorX_);
    unitVectorY_ /= mag(unitVectorY_);
    
    vector midPoint = startPoint_ + (startPoint_ - endPoint_)*0.5;
    
    startPointX_ = midPoint - dX_*unitVectorX_;
    endPointX_ = midPoint + dX_*unitVectorX_;
    
    startPointY_ = midPoint - dY_*unitVectorY_;
    endPointY_ = midPoint + dY_*unitVectorY_;
    
}



// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

twoDimBins::~twoDimBins()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

List<label> twoDimBins::isPointWithinBin
(
    const vector& rI,
    const label& cellI
)
{
    List<label> binNumbers;
binNumbers.append(-1);
binNumbers.append(-1);
    
    scalar rSEMag = mag(endPoint_ - startPoint_);

    vector rSI = rI - startPoint_;
    scalar rD = rSI & unitVector_;

    if
    (
        (rD >= 0) && (rD <= rSEMag)
    )
    {
        vector rSIX = rI - startPointX_;
        scalar rDX = rSIX & unitVectorX_;
        label nX = label(rDX/binWidthX_);

        vector rSIY = rI - startPointY_;
        scalar rDY = rSIY & unitVectorY_;
        label nY = label(rDY/binWidthY_);

        if( (nX >= 0) && (nX < nBinsX_) && (nY >= 0) && (nY < nBinsY_) )
        {
            binNumbers[0] = nX;
            binNumbers[1] = nY;
        }
    }


    return binNumbers;
}

scalarField twoDimBins::binPositionsX()
{
    scalarField positionsX(nBinsX_, 0.0);
    return positionsX;
}

scalarField twoDimBins::binPositionsY()
{
    scalarField positionsY(nBinsY_, 0.0);

    return positionsY;
}


void twoDimBins::write
(
    const fileName& path,
    const word& name
)
{

}

vector twoDimBins::position(/*const vector& h,*/ const scalar& r, const scalar& theta)
{
    vector p = vector::zero;
//     vector p = r*cos(theta)*angleUnitVectorY_ 
//                                     + r*sin(theta)*angleUnitVectorX_;
    return p;
}

const vector& twoDimBins::refVector() const
{
    return unitVector_;
}

List<label> twoDimBins::nBins()
{
    List<label> nBins;
nBins.append(-1);
nBins.append(-1);

//     nBins[0] = nBinsL_;
    nBins[0] = nBinsX_;
    nBins[1] = nBinsY_;

    return nBins;
}

scalar twoDimBins::binVolume(const label& n)
{
    return binWidthX_*binWidthY_*h_;
}


// * * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
