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
    collisionDistribution

Description

\*----------------------------------------------------------------------------*/

#include "collisionDistribution.H"
#include "graph.H"

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //
void collisionDistribution::setRadius()
{
    for(label i = 0; i < radius_.size(); i++)
    {
       radius_[i] = (0.5 + scalar(i)) * binWidth();
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Null constructor
collisionDistribution::collisionDistribution()
:
    distribution(),
    name_(),
    rMax_(-1),
    p_(),
    radius_()
{}

//- Construct from name
//- (used for reading an rdf from input file)
collisionDistribution::collisionDistribution(const word& name)
:
    distribution(),
    name_(name),
    rMax_(-1),
    p_(),
    radius_()
{}


// Construct from binWidth, rMax and name 
collisionDistribution::collisionDistribution
(
    const scalar& binWidth,
    const scalar& rMax,
    const word& name
)
:
    distribution(binWidth),
    name_(name),
    rMax_(rMax),
    p_(label((rMax/binWidth) + 1), 0.0),
    radius_(label((rMax/binWidth) + 1), 0.0)
{
    setRadius();
}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

collisionDistribution::~collisionDistribution()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void collisionDistribution::addPoint(const scalar& r)
{
    if(r < rMax_)
    {
        add(r);
    }
}

void collisionDistribution::setDistribution()
{
    //-reset probabilities
    forAll(p_, bin)
    {
        p_[bin] = 0.0;
    }

    iterator iter(this->begin());

    //- make sure that we have a quantity at r = 0
    scalar valueToAdd = 0.0;

    label n = label(valueToAdd/binWidth()) - label(neg(valueToAdd/binWidth()));

    iter = find(n);

    if (iter == this->end())
    {
        this->insert(n,0);
    }

//     List< Pair<scalar> > data = raw(); //normalised()

    insertMissingKeys();

    List<label> keys = toc();

    sort(keys);

    List< Pair<scalar> > normDist(size());

    scalar maxValue = 0.0;

    forAll(keys, k)
    {
        if(k < p_.size())
        {
            label key = keys[k];
    
            scalar nMols = scalar((*this)[key]);
    
            if( nMols > maxValue)
            {
                maxValue = nMols;
            }
        }
    }

    forAll(keys, k)
    {
        if(k < p_.size())
        {
            label key = keys[k];

            scalar scaledValue = 0.0;

            if(maxValue > 0.0)
            {
                scaledValue = 1/maxValue;
            }

            p_[k] = scalar((*this)[key])*scaledValue;
        }
    }
/*
    forAll(data, k)
    {
        if(k < p_.size())
        {
            p_[k] = data[k].second();
        }
    }*/

//     clear();
}

// bin value for g
// scalar collisionDistribution::gBin(const scalar& r) const
// {
//     label key = label(r/binWidth());
// 
//     if(r < 0.0)
//     {
//         return 0.0;
//     }
//     else if
//     (
//         (key < g_.size()) &&
//         (key >= 0)
//     )
//     {
//         return g_[key];
//     }
//     else
//     {
//         return 1.0;
//     }
// }

// linear interpolation for g
// scalar collisionDistribution::gLinear(const scalar& r) const
// {
//     label key = label(r/binWidth());
// 
//     if(r < 0.0)
//     {
//         return 0.0;
//     }
//     else if
//     (
//         (key < g_.size()) &&
//         (key >= 0)
//     )
//     {
//         scalar deltaR = radius_[key] - r;
// 
//         label key2 = key;
// 
//         if(deltaR < 0.0)
//         {
//             key2++;
//         }
//         else
//         {
//             key2--;
//         }
// 
//         const scalar& g1 = g_[key];
//         const scalar& g2 = g_[key2];
// 
//         scalar gNew = ( ((g2 - g2)/binWidth()) * (mag(deltaR)) ) + g1;
// 
//         return gNew;
//     }
//     else
//     {
//         return 0.0;
//     }
// }


//- return col. distr. in the form of a list of pairs (p and r)
List< Pair<scalar> > collisionDistribution::distrib()
{
    List< Pair<scalar> > radialDist(p_.size());

    forAll(p_, bin)
    {
        radialDist[bin].first() = radius_[bin];

        radialDist[bin].second() = p_[bin];
    }

    return radialDist;
}



void collisionDistribution::clearRadialDistribution()
{
    // clear the list
    forAll(p_, bin)
    {
        p_[bin] = 0.0;
    }

    //clear the hash table
    clear();
}

void collisionDistribution::writeDistribution
(
    const Time& runTime
)
{
    if(runTime.outputTime())
    {
        fileName timePath(runTime.path()/runTime.timeName()/"uniform");

        if (!isDir(timePath))
        {
            mkDir(timePath);
        }

        OFstream colDistrFile(timePath/name_+".raw");

        if (colDistrFile.good())
        {
            colDistrFile << distrib() << endl;
        }
        else
        {
            FatalErrorIn("void collisionDistribution::writeRDF()")
                << "Cannot open file " << colDistrFile.name()
                << abort(FatalError);
        }

        writeTimeData(timePath, name_, radius_ , p_);
    }
}

// void collisionDistribution::writeTimeData
// (
//     const fileName& pathName,
//     const word& nameFile,
//     const scalarField& xData,
//     const scalarField& yData
// )
// {
//     fileName writeFile(pathName/nameFile);
// 
//     graph outputGraph("title", "x", "y", xData, yData);
// 
//     outputGraph.write(writeFile, "raw");
// }


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

// void collisionDistribution::operator=(const collisionDistribution& rhs)
// {
//     // Check for assignment to self
//     if (this == &rhs)
//     {
//         FatalErrorIn("collisionDistribution::operator=(const collisionDistribution&)")
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

// Ostream& operator<<(Ostream& os, const collisionDistribution& d)
// {
//     os  << d.binWidth_
//         << static_cast<const Map<label>&>(d);
// 
//     // Check state of Ostream
//     os.check
//     (
//         "Ostream& operator<<(Ostream&, "
//         "const collisionDistribution&)"
//     );
// 
//     return os;
// }

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
