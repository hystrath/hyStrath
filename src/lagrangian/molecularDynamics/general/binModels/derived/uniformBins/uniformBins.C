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
    uniformBins

Description

\*----------------------------------------------------------------------------*/

#include "uniformBins.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(uniformBins, 0);

addToRunTimeSelectionTable(binModel, uniformBins, dictionary);



// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

//- Construct from components
uniformBins::uniformBins
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
    rSEMag_(mag(endPoint_ - startPoint_)),
    nBins_(readLabel(propsDict_.lookup("nBins"))),
    binWidth_(mag(endPoint_ - startPoint_)/(nBins_)),
    area_(readScalar(propsDict_.lookup("area")))
{
//     Info << "rSEMag = " << rSEMag_ << endl;
//     Info << "startPoint = " << startPoint_ << endl;
//     Info << "endPoint = " << endPoint_ << endl;
}



// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

uniformBins::~uniformBins()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// cellI is a dummy variable
label uniformBins::isPointWithinBin
(
    const vector& rI,
    const label& cellI
)
{
    label binNumber = -1;

//     scalar rSEMag = mag(endPoint_ - startPoint_);

    vector rSI = rI - startPoint_;
    scalar rD = rSI & unitVector_;
    label n = label(rD/binWidth_);

//     Info << "n = " << n << ", rD = " << rD 
//          << endl;
    
    
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


scalarField uniformBins::binPositions()
{
    scalarField positions(nBins_, 0.0);

    forAll(positions, i)
    {
        positions[i] = 0.5*binWidth_ + scalar(i)*binWidth_;
    }

    return positions;
}

vectorField uniformBins::bins()
{
    vectorField positions(nBins_, vector::zero);

    forAll(positions, i)
    {
        positions[i] = startPoint_ + (0.5 + scalar(i))*binWidth_*unitVector_;
    }

    return positions;
}

const label& uniformBins::nBins() const
{
    return nBins_;
}

scalar uniformBins::binVolume(const label& n)
{
    scalar volume = area_*binWidth_;

    return volume;
}


// * * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
