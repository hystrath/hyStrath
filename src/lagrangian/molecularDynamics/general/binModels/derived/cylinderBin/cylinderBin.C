/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016-2020 hyStrath
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of hyStrath, a derivative work of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    cylinderBin

Description

\*----------------------------------------------------------------------------*/

#include "cylinderBin.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(cylinderBin, 0);

addToRunTimeSelectionTable(binModel, cylinderBin, dictionary);



// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

//- Construct from components
cylinderBin::cylinderBin
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
    rMaxSqr_(sqr(readScalar(propsDict_.lookup("radius")))),
    area_(0.0)
{
    area_ = rMaxSqr_*constant::mathematical::pi;
}



// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

cylinderBin::~cylinderBin()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

label cylinderBin::isPointWithinBin
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
            vector rBin = startPoint_ + scalar(n)*binWidth_*unitVector_;
            vector rBI = rBin - rI;

            scalar distSqr = magSqr(rBI) - magSqr(rBI & unitVector_);

            if(distSqr <= rMaxSqr_)
            {
                binNumber = n;
            }
        }
    }

    return binNumber;
}


scalarField cylinderBin::binPositions()
{
    scalarField positions(nBins_, 0.0);

    forAll(positions, i)
    {
        positions[i] = 0.5*binWidth_ + scalar(i)*binWidth_;
    }

    return positions;
}

vectorField cylinderBin::bins()
{
    vectorField positions(nBins_, vector::zero);

    forAll(positions, i)
    {
        positions[i] = startPoint_ + (0.5 + scalar(i))*binWidth_*unitVector_;
    }

    return positions;
}

const label& cylinderBin::nBins() const
{
    return nBins_;
}

scalar cylinderBin::binVolume(const label& n)
{
    scalar volume = area_*binWidth_;

    return volume;
}


// * * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
