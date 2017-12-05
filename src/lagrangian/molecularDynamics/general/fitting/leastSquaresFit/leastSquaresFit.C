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
    leastSquaresFit

Description

\*----------------------------------------------------------------------------*/

#include "leastSquaresFit.H"
#include "graph.H"

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //
// void leastSquaresFit::setRadius()
// {
//     for(label i = 0; i < noOfBins_; i++)
//     {
//        radius_[i] = (0.5 + scalar(i)) * binWidth();
//     }
// }

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Null constructor
leastSquaresFit::leastSquaresFit()
:
    name_(),
    write_(),
    xs_(),
    ys_(),
    gradient_(0.0),
    yIntercept_(0.0)
{}



leastSquaresFit::leastSquaresFit
(
    const word& name,
    const bool& write,
    const scalarField& xs,
    const scalarField& ys
)
:
    name_(name),
    write_(write),
    xs_(xs),
    ys_(ys),
    gradient_(0.0),
    yIntercept_(0.0)
{

    if(xs_.size() != ys_.size())
    {
        FatalErrorIn("leastSquaresFit")
            << "size of x and y lists are inconsistent"
            << abort(FatalError);
    }

    setFitParameters();
}

void leastSquaresFit::setFitParameters()
{
    scalar N = scalar(xs_.size());


    scalar termA = 0.0;
    scalar termB = 0.0;
    scalar termC = 0.0;
    scalar termD = 0.0;

    forAll(xs_, n)
    {
        termA += xs_[n]*xs_[n];
        termB += xs_[n];
        termC += xs_[n]*ys_[n];
        termD += ys_[n];
    }

    //solve simultaneous equations

    yIntercept_ = (termC - ((termA*termD)/termB)) /
                    (termB-((termA*N)/termB));

    gradient_ = (termD- (N*yIntercept_) )/termB;


    if(yIntercept_ > 0)
    {
        Info << "y = " << gradient_ <<"x + " << yIntercept_ << endl;
    }
    else if (yIntercept_ < 0)
    {
        Info << "y = " << gradient_ <<"x - " << mag(yIntercept_) << endl;
    }
    else
    {
        Info << "y = " << gradient_ <<"x" << endl;
    }
}

void leastSquaresFit::setInitialData
(
    const word& name,
    const bool& write,
    const label nBins
)
{
    name_ = name;
    write_ = write;
    xs_.setSize(nBins, 0.0);
    ys_.setSize(nBins, 0.0);

//     Info << "nBins: " << nBins << endl;
}

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

leastSquaresFit::~leastSquaresFit()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //



void leastSquaresFit::writeField
(
    const Time& runTime
)
{
    if(runTime.outputTime())
    {
        if(write_)
        {
            fileName timePath(runTime.path()/runTime.timeName()/"uniform");
    
            if (!isDir(timePath))
            {
                mkDir(timePath);
            }
    
            writeTimeData(timePath, name_, xs_ , ys_);
        }
    }
}

void leastSquaresFit::writeTimeData
(
    const fileName& pathName,
    const word& nameFile,
    const scalarField& xData,
    const scalarField& yData
)
{
    fileName writeFile(pathName/nameFile);

    graph outputGraph("title", "x", "y", xData, yData);

    outputGraph.write(writeFile, "raw");
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

// void leastSquaresFit::operator=(const leastSquaresFit& rhs)
// {
//     // Check for assignment to self
//     if (this == &rhs)
//     {
//         FatalErrorIn("leastSquaresFit::operator=(const leastSquaresFit&)")
//             << "Attempted assignment to self"
//             << abort(FatalError);
//     }
// 
//     Map<label>::operator=(rhs);
// 
//     binWidth_ = rhs.binWidth();
// }


// * * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

// Ostream& operator<<(Ostream& os, const leastSquaresFit& d)
// {
//     os  << d.binWidth_
//         << static_cast<const Map<label>&>(d);
// 
//     // Check state of Ostream
//     os.check
//     (
//         "Ostream& operator<<(Ostream&, "
//         "const leastSquaresFit&)"
//     );
// 
//     return os;
// }

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
