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
    simpleBins

Description

\*----------------------------------------------------------------------------*/

#include "simpleBins.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(simpleBins, 0);

addToRunTimeSelectionTable(binModel, simpleBins, dictionary);



// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //



//- Construct from components
simpleBins::simpleBins
(
    const polyMesh& mesh,
    const dictionary& dict
)
:
    binModel(mesh, dict),
    propsDict_(dict.subDict(typeName + "Properties")),
    bb_(mesh.bounds().min(), mesh.bounds().max()),
    nBins_(readLabel(propsDict_.lookup("nBins")))  
{
    
    const word option = propsDict_.lookup("direction");
    
    if(option == "x")
    {
        unitVector_ = vector(1,0,0);
        area_ = (bb_.span() & vector(0,1,0)) * (bb_.span() & vector(0,0,1));
    }
    else if(option == "y")
    {
        unitVector_ = vector(0,1,0);
      
        area_ = (bb_.span() & vector(1,0,0)) * (bb_.span() & vector(0,0,1));
        
    }
    else if(option == "z")
    {
        unitVector_ = vector(0,0,1);
       
        area_ = (bb_.span() & vector(0,1,0)) * (bb_.span() & vector(1,0,0));
        
    }
    else
    {
        // Error
        
    }
    
    startPoint_ = bb_.min();
    endPoint_ = (bb_.span() & unitVector_)*unitVector_ + bb_.min();
    rSEMag_ = mag(endPoint_ - startPoint_);
    binWidth_ = rSEMag_/scalar(nBins_);
    
    Info << "start point = " << startPoint_ << endl;
    Info << "end point = " << endPoint_ << endl;
    Info << "rSEMag = " << rSEMag_ << endl;    
    
}



// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

simpleBins::~simpleBins()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// cellI is a dummy variable
label simpleBins::isPointWithinBin
(
    const vector& rI,
    const label& cellI
)
{
    label binNumber = -1;

    vector rSI = rI - startPoint_;
    scalar rD = rSI & unitVector_;
    label n = label(rD/binWidth_);

//     Info << "n = " << n << endl;
    if
    (
        (rD <= rSEMag_) && (rD >= 0.0)
    )
    {

        if(n == nBins_)
        {
            n--;
        }

        if( (n >=0) && (n < nBins_) )
        {
            binNumber = n;
        }
    }

    return binNumber;
}


scalarField simpleBins::binPositions()
{
    scalarField positions(nBins_, 0.0);

    forAll(positions, i)
    {
        positions[i] = 0.5*binWidth_ + scalar(i)*binWidth_;
    }

    return positions;
}

vectorField simpleBins::bins()
{
    vectorField positions(nBins_, vector::zero);

    forAll(positions, i)
    {
        positions[i] = startPoint_ + (0.5 + scalar(i))*binWidth_*unitVector_;
    }

    return positions;
}

const label& simpleBins::nBins() const
{
    return nBins_;
}

scalar simpleBins::binVolume(const label& n)
{
    scalar volume = area_*binWidth_;

    return volume;
}


// * * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
