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
    sphericalBins

Description

\*----------------------------------------------------------------------------*/

#include "sphericalBins.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(sphericalBins, 0);

addToRunTimeSelectionTable(binModel, sphericalBins, dictionary);



// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

label sphericalBins::findBin(const scalar& r)
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

        scalar rLimit1 = magRadii_[n] - 0.5*binWidths_[n];
        scalar rLimit2 = magRadii_[n] + 0.5*binWidths_[n];
    
        if((r >= rLimit1) && (r < rLimit2))
        {}
        else
        {
            n++;
        }
    }

    return n;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

//- Construct from components
sphericalBins::sphericalBins
(
    const polyMesh& mesh,
    const dictionary& dict
)
:
    binModel(mesh, dict),
    propsDict_(dict.subDict(typeName + "Properties")),
    startPoint_(propsDict_.lookup("startPoint")),
    radius_(readScalar(propsDict_.lookup("radius"))),
    nBins_(readLabel(propsDict_.lookup("nBins"))),
    binWidth_(radius_/scalar(nBins_)),

    magRadii_(),
    binWidths_(),
    avVolume_(0.0),
    minBinWidth_(0.0),
    n_()
{
    //- set volumes --- here we impose a constant volume in all bins,
    // hence the binWidth is going change radially

    avVolume_ = 4.0*radius_*radius_*radius_*constant::mathematical::pi/(3.0*nBins_);

    Info << "average volume: " << avVolume_ << endl;

    DynamicList<scalar> radii(0);
    DynamicList<scalar> binWidths(0);

    scalar binWidth = Foam::pow((3.0*avVolume_/(4.0*constant::mathematical::pi)), (1.0/3.0));

    Info << "starting binWidth: " << binWidth << endl;

    binWidths.append(binWidth);

    radii.append(0.5*binWidths[0]);

    for(label i = 1; i < nBins_; i++)
    {
        scalar r2 = radii[i-1] + 0.5*binWidths[i-1];
        scalar r1 = pow(((3.0*avVolume_/(4.0*constant::mathematical::pi)) + (r2*r2*r2)), (1.0/3.0));

        if(r2 <= radius_)
        {
            radii.append(0.5*(r1+r2));
            binWidths.append(r1-r2);
        }
        else
        {
            break;
        }
    }

    //binWidths.shrink();
    //radii.shrink();

    nBins_ = binWidths.size();

    magRadii_.setSize(nBins_, 0.0);
    binWidths_.setSize(nBins_, 0.0);

    forAll(magRadii_, n)
    {
        magRadii_[n] = radii[n];
        binWidths_[n] = binWidths[n];
    }

    minBinWidth_ = binWidths_[nBins_-1]/3.0;

    n_.setSize(label(radius_/minBinWidth_), 0);

    forAll(n_, n)
    {
        scalar r = 0.5*minBinWidth_ + scalar(n)*minBinWidth_;

        for(label i = 0; i < nBins_; i++)
        {
            scalar rLimit1 = magRadii_[i] - 0.5*binWidths_[i];
            scalar rLimit2 = magRadii_[i] + 0.5*binWidths_[i];

            if((r >= rLimit1) && (r < rLimit2))
            {
                n_[n] = i;
            }
        }
    }

    Info << "binWidths_: " << binWidths_ << endl;
}



// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

sphericalBins::~sphericalBins()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //






label sphericalBins::isPointWithinBin
(
    const vector& rI,
    const label& cellI
)
{
    label binNumber = -1;

    scalar rSIMag = mag(rI - startPoint_);

    if(rSIMag <= radius_)
    {
        label n = findBin(rSIMag);
    
        if
        (
            n != -1
        )
        {
            binNumber = n;
        }
    }

    return binNumber;
}


scalarField sphericalBins::binPositions()
{
//     scalarField positions(nBins_, 0.0);
// 
//     forAll(positions, i)
//     {
//         positions[i] = 0.5*binWidth_ + scalar(i)*binWidth_;
//     }

    return magRadii_;
}


// this is a work around
vectorField sphericalBins::bins()
{
    vectorField positions(nBins_, vector::zero);

    forAll(positions, i)
    {
        positions[i] = magRadii_[i]*vector(1,0,0);
    }

    return positions;
}

const label& sphericalBins::nBins() const
{
    return nBins_;
}

scalar sphericalBins::binVolume(const label& n)
{
    return avVolume_;
}


// * * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
