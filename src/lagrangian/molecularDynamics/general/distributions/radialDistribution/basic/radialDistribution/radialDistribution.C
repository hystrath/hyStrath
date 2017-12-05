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
    radialDistribution

Description

\*----------------------------------------------------------------------------*/

#include "radialDistribution.H"
#include "graph.H"

// using namespace Foam::constant::mathematical;

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //
void radialDistribution::setRadius()
{
    for(label i = 0; i < noOfBins_; i++)
    {
       radius_[i] = (0.5 + scalar(i)) * binWidth();
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Null constructor
radialDistribution::radialDistribution()
:
//     distribution(),
    name_(),
    fileName_("radialDistributions"),
    rMax_(-1),
    noOfBins_(-1),
    binWidth_(-1.0),
    g_(),
    mols_(),
    radius_()
{}

//- Construct from name
//- (used for reading an rdf from input file)
radialDistribution::radialDistribution(const word& name)
:
//     distribution(),
    name_(name),
    fileName_("radialDistributions"),
    rMax_(-1.0),
    noOfBins_(-1),
    binWidth_(-1.0),
    g_(),
    mols_(),
    radius_()
{}


//- Construct from name, rMax and noOfBins 
//- (used for sampling radial distribution)
radialDistribution::radialDistribution
(
    const word& name,
    const scalar& rMax,
    const label& noOfBins
)
:
//     distribution(binWidth),
    name_(name),
    fileName_("radialDistributions"),
    rMax_(rMax),
    noOfBins_(noOfBins),
    binWidth_(rMax/noOfBins_),
    g_(noOfBins, 0.0),
    mols_(noOfBins, 0.0),
    radius_(noOfBins, 0.0)
{
    setRadius();
}

// Construct from components
//- initial conditions or for cases in which rdf remains fixed
radialDistribution::radialDistribution
(
    const word& name,
    const List< scalar >& g,
    const List< scalar >& radius
)
:
//     distribution(binWidth),
    name_(name),
    fileName_("radialDistributions"),
    rMax_(radius[g.size()-1]),
    noOfBins_(g.size()),
    binWidth_(-1.0),
    g_(g),
    mols_(noOfBins_, 0.0),
    radius_(radius)
{
    binWidth_ = readBinWidth();
}

//- construct from components: pair lists
radialDistribution::radialDistribution
(
    const word& name,
    const List< Pair<scalar> >& rdf
//     const scalar& binWidth
)
:
    name_(name),
    fileName_("radialDistributions"),
    rMax_(-1.0),
    noOfBins_(rdf.size()),
    binWidth_(-1.0),
    g_(),
    mols_(noOfBins_, 0.0),
    radius_()

{
    setRDF(rdf);
}


void radialDistribution::setRdfName
(
    const word& name
)
{
    name_ = name;
}

void radialDistribution::setRDF(const List< Pair<scalar> >& rdf)
{
    g_.setSize(rdf.size());
    radius_.setSize(rdf.size());

    forAll(rdf, bin)
    {
        radius_[bin] = rdf[bin].first();

        g_[bin] = rdf[bin].second();
    }

    binWidth_ = readBinWidth();

    rMax_ = rdf[rdf.size()-1].first() + 0.5*binWidth_;
}

void radialDistribution::setRdf
(
    const word& name,
    const List< scalar >& g,
    const List< scalar >& radius
)
{
    name_ = name;

    g_.setSize(g.size());
    g_ = g;

    radius_.setSize(radius.size());
    radius_ = radius;

    binWidth_ = readBinWidth();

    noOfBins_ = g.size();
    rMax_ = radius[noOfBins_-1] + 0.5*binWidth_;
}



void radialDistribution::setRDF
(
//     const scalar& binWidth,
    const scalar& rMax,
    const label& noOfBins
)
{
    rMax_ = rMax;
    noOfBins_ = noOfBins;

    binWidth_ = rMax/noOfBins_,

    g_.setSize(noOfBins, 0.0);

    radius_.setSize(noOfBins, 0.0);

    setRadius();
}


scalar radialDistribution::readBinWidth()
{
    //- test for zero size
    if(g_.size() == 0)
    {
        FatalErrorIn("radialDistribution")
            << "Size of RDF:  " << name_
            << " is zero" 
            << abort(FatalError);
    }

    //- test for common bin width

    scalar binWidth = 0.0;

    bool constantBinWidth = true;

    for(label i = 0; i < radius_.size() - 1; i++)
    {
        if(binWidth == 0.0)
        {
            binWidth = mag(radius_[i+1] - radius_[i]);
//             Info << " width: " << binWidth << endl;
        }
        else
        {
            scalar newBinWid = mag(radius_[i+1] - radius_[i]);

//             Info << " width: " << newBinWid << endl;

            if(mag(newBinWid - binWidth) > SMALL) 
            {
                constantBinWidth = false;
                break;
            }
        }
    }

    if(!(binWidth > SMALL) || !constantBinWidth)
    {
        FatalErrorIn("radialDistribution")
            << "Check binWidth in RDF: " << name_
            << " for constant binWidth" 
            << abort(FatalError);
    }

    return binWidth;
}
// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

radialDistribution::~radialDistribution()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void radialDistribution::addToDistribution(const scalar& r)
{
    label n = label(r/binWidth_);

    if
    (
        (n < g_.size()) &&
        (n >= 0)
    )
    {
        mols_[n] += 1.0;
    }
}

void radialDistribution::setRdf(const scalar& rhoAv, const scalar& nMols)
{

    forAll(g_, r)
    {
        const scalar& radius = radius_[r];

        //- volume of a shell
        scalar Vshell = 4*constant::mathematical::pi*binWidth_*sqr(radius);

        scalar scaledValue = 1/(Vshell*rhoAv*nMols);

        g_[r] = mols_[r]*scaledValue;
    }
}

// bin value for g
scalar radialDistribution::gBin(const scalar& r) const
{
    label key = label(r/binWidth_);

    if
    (
        (key < g_.size()) &&
        (key >= 0)
    )
    {
        return g_[key];
    }
    else
    {
        return 0.0;
    }
}

// linear interpolation for g
scalar radialDistribution::gLinear(const scalar& r) const
{
    label key = label(r/binWidth());

    if
    (
        (key < g_.size()) &&
        (key >= 0)
    )
    {
        scalar deltaR = radius_[key] - r;

        label key2 = key;

        if(deltaR < 0.0)
        {
            key2++;
        }
        else
        {
            key2--;

            if(key2 < 0)
            {
                key2 = 0;
            }
        }

        const scalar& g1 = g_[key];
        const scalar& g2 = g_[key2];

        scalar gNew = ( ((g2 - g2)/binWidth()) * (mag(deltaR)) ) + g1;

        return gNew;
    }
    else
    {
        return 0.0;
    }
}


//- return rdf in the form of a list of pairs (g and r)
List< Pair<scalar> > radialDistribution::radialDistrib()
{
    List< Pair<scalar> > radialDist(g_.size());

    forAll(g_, bin)
    {
        radialDist[bin].first() = radius_[bin];

        radialDist[bin].second() = g_[bin];
    }

    return radialDist;
}

void radialDistribution::clearCollector()
{
    // clear the list
    mols_ = 0.0;
}

void radialDistribution::clearRadialDistribution()
{
    // clear the list
    g_ = 0.0;
    clearCollector();

    //clear the hash table
//     clear();
}

void radialDistribution::writeRDF
(
    const Time& runTime
)
{
    if(runTime.outputTime())
    {
        fileName timePath(runTime.path()/runTime.timeName()/"uniform"/fileName_);

        if (!isDir(timePath))
        {
            mkDir(timePath);
        }

        OFstream rdfFile(timePath/name_+".raw");

        if (rdfFile.good())
        {
            rdfFile << radialDistrib() << endl;
        }
        else
        {
            FatalErrorIn("void radialDistribution::writeRDF()")
                << "Cannot open file " << rdfFile.name()
                << abort(FatalError);
        }

        writeTimeData(timePath, name_, radius_ , g_);
//         clearRDF();
    }
}

void radialDistribution::writeTimeData
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

void radialDistribution::readRDF
(
    const Time& runTime
)
{
    fileName timePath(runTime.path()/runTime.timeName()/"uniform"/fileName_);

    IFstream rdfFile(timePath/name_+".raw");

    List< Pair<scalar> > rdf;

    if (rdfFile.good())
    {
        rdfFile >> rdf;
    }
    else
    {
        FatalErrorIn("void radialDistribution::readRDF()")
            << "Cannot open file " << rdfFile.name()
            << abort(FatalError);
    }

    setRDF(rdf);

    // set bin width of distribution class (this)
//     setBinWidth(readBinWidth()); 
}

// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

// void radialDistribution::operator=(const radialDistribution& rhs)
// {
//     // Check for assignment to self
//     if (this == &rhs)
//     {
//         FatalErrorIn("radialDistribution::operator=(const radialDistribution&)")
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

// Ostream& operator<<(Ostream& os, const radialDistribution& d)
// {
//     os  << d.binWidth_
//         << static_cast<const Map<label>&>(d);
// 
//     // Check state of Ostream
//     os.check
//     (
//         "Ostream& operator<<(Ostream&, "
//         "const radialDistribution&)"
//     );
// 
//     return os;
// }

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
