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
    cuboidBin

Description

\*----------------------------------------------------------------------------*/

#include "cuboidBin.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(cuboidBin, 0);

addToRunTimeSelectionTable(binModel, cuboidBin, dictionary);



// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

//- Construct from components
cuboidBin::cuboidBin
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
}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

cuboidBin::~cuboidBin()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


label cuboidBin::isPointWithinBin
(
    const vector& rI,
    const label& cellI
)
{
    label binNumber = -1;

    scalar rSEMag = mag(endPoint_ - startPoint_);

    vector rSI = rI - startPoint_;
    scalar rD = rSI & unitVector_;
    label n = label(rD/binWidth_);

    if
    (
        (n >= 0) && (rD <= rSEMag)
    )
    {

        if(n == nBins_) 
        {
            n--;
        }

        if(n < nBins_)
        {
            vector rBin = startPoint_ + (0.5+scalar(n))*binWidth_*unitVector_;
            vector rBI = rBin - rI;
    
            scalar dist1 = mag(rBI & n1_);
    
            if(dist1 <= d1_)
            {
                scalar dist2 = mag(rBI & n2_);
    
                if(dist2 <= d2_)
                {
                    binNumber = n;
                }
            }
        }
    }

    return binNumber;
}


scalarField cuboidBin::binPositions()
{
    scalarField positions(nBins_, 0.0);

    forAll(positions, i)
    {
        positions[i] = 0.5*binWidth_ + scalar(i)*binWidth_;
    }

    return positions;
}

vectorField cuboidBin::bins()
{
    vectorField positions(nBins_, vector::zero);

    forAll(positions, i)
    {
        positions[i] = startPoint_ + (0.5 + scalar(i))*binWidth_*unitVector_;
    }

    return positions;
}

const label& cuboidBin::nBins() const
{
    return nBins_;
}

scalar cuboidBin::binVolume(const label& n)
{
    scalar volume = area_*binWidth_;

    return volume;
}


// * * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
