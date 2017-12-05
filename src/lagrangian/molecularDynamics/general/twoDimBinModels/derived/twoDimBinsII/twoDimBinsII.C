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
    twoDimBinsII

Description

\*----------------------------------------------------------------------------*/

#include "twoDimBinsII.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(twoDimBinsII, 0);

addToRunTimeSelectionTable(twoDimBinModel, twoDimBinsII, dictionary);



// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

//- Construct from components
twoDimBinsII::twoDimBinsII
(
    const polyMesh& mesh,
    const dictionary& dict
)
:
    twoDimBinModel(mesh, dict),
    propsDict_(dict.subDict(typeName + "Properties")),
    startPoint_(propsDict_.lookup("startPoint")),
    endPoint_(propsDict_.lookup("endPoint")),
    unitVectorX_(propsDict_.lookup("unitVectorX")),
    unitVectorY_(propsDict_.lookup("unitVectorY")),
    unitVectorZ_(propsDict_.lookup("unitVectorZ")),
    nBinsX_(readLabel(propsDict_.lookup("nBinsX"))),
    nBinsY_(readLabel(propsDict_.lookup("nBinsY")))
{

    unitVectorX_ /= mag(unitVectorX_);
    unitVectorY_ /= mag(unitVectorY_);
    unitVectorZ_ /= mag(unitVectorZ_);  
    box_.resetBoundedBox(startPoint_, endPoint_);
    
    vector rS = box_.span();
    
    lengthX_ = rS & unitVectorX_;
    lengthY_ = rS & unitVectorY_;
    lengthZ_ = rS & unitVectorZ_; 
    
    binWidthX_ = lengthX_/nBinsX_;
    binWidthY_ = lengthY_/nBinsY_;
    
//     unitVectorZ_ = unitVectorX_ ^ unitVectorY_;
//     unitVectorZ_ /= mag(unitVectorZ_);  
  
 
    
    Info  << nl << "twoDimBinsII properties" << nl 
        << "binWidthX: " << binWidthX_ << nl
        << "binWidthY: " << binWidthY_
        << endl;
    
}



// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

twoDimBinsII::~twoDimBinsII()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

List<label> twoDimBinsII::isPointWithinBin
(
    const vector& rI,
    const label& cellI
)
{
    List<label> binNumbers;
binNumbers.append(-1);
binNumbers.append(-1);
    
    if(box_.contains(rI))
    {
        vector rSI = rI - startPoint_;
        scalar rDx = rSI & unitVectorX_;
        scalar rDy = rSI & unitVectorY_;
        label nX = label(rDx/binWidthX_);
        label nY = label(rDy/binWidthY_);
        
        if( (nX >= 0) && (nY >= 0) )
        {
            if(nX == nBinsX_) 
            {
                nX--;
            }
            
            binNumbers[0] = nX;
        
            if(nY == nBinsY_) 
            {
                nY--;
            }
            
            binNumbers[1] = nY;            
        }    
    }
    
    return binNumbers;
}

scalarField twoDimBinsII::binPositionsX()
{
    scalarField positionsX(nBinsX_, 0.0);
    
    forAll(positionsX, i)
    {
        positionsX[i] = binWidthX_*i + binWidthX_*0.5;
    }
    
    return positionsX;
}

scalarField twoDimBinsII::binPositionsY()
{
    scalarField positionsY(nBinsY_, 0.0);

    forAll(positionsY, i)
    {
        positionsY[i] = binWidthY_*i + binWidthY_*0.5;
    }    
    
    return positionsY;
}


void twoDimBinsII::write
(
    const fileName& path,
    const word& name
)
{

}

vector twoDimBinsII::position(/*const vector& h,*/ const scalar& r, const scalar& theta)
{
    vector p = vector::zero;
//     vector p = r*cos(theta)*angleUnitVectorY_ 
//                                     + r*sin(theta)*angleUnitVectorX_;
    return p;
}

const vector& twoDimBinsII::refVector() const
{
    return unitVectorX_;
}

List<label> twoDimBinsII::nBins()
{
    List<label> nBins;
nBins.append(-1);
nBins.append(-1);

//     nBins[0] = nBinsL_;
    nBins[0] = nBinsX_;
    nBins[1] = nBinsY_;

    return nBins;
}

scalar twoDimBinsII::binVolume(const label& n)
{
    return binWidthX_*binWidthY_*lengthZ_;
}


// * * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
