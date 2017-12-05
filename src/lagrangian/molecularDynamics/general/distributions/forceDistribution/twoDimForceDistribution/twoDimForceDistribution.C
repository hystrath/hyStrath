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
    twoDimForceDistribution

Description

\*----------------------------------------------------------------------------*/

#include "twoDimForceDistribution.H"
#include "graph.H"
#include "IFstream.H"
#include "OFstream.H"

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //
// void twoDimForceDistribution::setRadii()
// {
//     for(label i = 0; i < noOfBins_; i++)
//     {
//         magRadii_[i] = 0.5*binWidth_ + scalar(i)*binWidth_;
// 
//         radii_[i] = startPoint_ + (0.5*binWidth_ + scalar(i)*binWidth_)*unitVector_;
//     }
// }

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Null constructor
twoDimForceDistribution::twoDimForceDistribution()
:
    name_(),
    fileName_("forceDistributions"),
    startPoint_(vector::zero),
    lengthX_(-1.0),
    lengthY_(-1.0),
    unitVectorX_(vector::zero),
    unitVectorY_(vector::zero),
    nBinsX_(-1),
    nBinsY_(-1),
    binWidthX_(0.0),
    binWidthY_(0.0),
    mols_(),
    forces_(),
    energies_()
{}





// Construct from dictionary
twoDimForceDistribution::twoDimForceDistribution
(
    const dictionary& dict
)
:
    name_(dict.lookup("distributionName")),
    fileName_("forceDistributions"),
    startPoint_(dict.lookup("startPoint")),
    lengthX_(readScalar(dict.lookup("lengthX"))),
    lengthY_(readScalar(dict.lookup("lengthY"))),
    unitVectorX_(dict.lookup("unitVectorX")),
    unitVectorY_(dict.lookup("unitVectorY")),
    nBinsX_(readLabel(dict.lookup("nBinsX"))),
    nBinsY_(readLabel(dict.lookup("nBinsY"))),
    binWidthX_(0.0),
    binWidthY_(0.0),
    mols_(nBinsX_),
    forces_(nBinsX_),
    energies_(nBinsX_)
{
    unitVectorX_ /= mag(unitVectorX_);
    unitVectorY_ /= mag(unitVectorY_);

    binWidthX_ = lengthX_/nBinsX_;
    binWidthY_ = lengthY_/nBinsY_;

    forAll(mols_, x)
    {
        mols_[x].setSize(nBinsY_, 0.0);
        forces_[x].setSize(nBinsY_, vector::zero);
        energies_[x].setSize(nBinsY_, 0.0);
    }

//     setRadii();
}


// Construct as copy
// twoDimForceDistribution::twoDimForceDistribution(const twoDimForceDistribution& d)
// :
//     Map<label>(static_cast< Map<label> >(d)),
//     binWidth_(d.binWidth())
// {}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

twoDimForceDistribution::~twoDimForceDistribution()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void twoDimForceDistribution::addToDistribution
(
    const vector& r, 
    const vector& force,
    const scalar& energy
)
{
    vector rLocal = r - startPoint_;

    scalar rX = rLocal & unitVectorX_;
    scalar rY = rLocal & unitVectorY_;

    if
    (
         (binWidthX_ > 0.0) &&
         (binWidthY_ > 0.0)
    )
    {
        label nX = label(rX/binWidthX_);
        label nY = label(rY/binWidthY_);
    
        if
        (
            (nX < nBinsX_) &&
            (nY < nBinsY_) &&
            (nX >= 0) &&
            (nY >= 0)
        )
        {
            forces_[nX][nY] += force;
            energies_[nX][nY] += energy;
            mols_[nX][nY] += 1.0;
        }
    }
}

bool twoDimForceDistribution::isWithinDistributionRange(const vector& r)
{
    vector rLocal = r - startPoint_;

    scalar rX = rLocal & unitVectorX_;
    scalar rY = rLocal & unitVectorY_;

    bool inRange = false;

    if
    (
         (binWidthX_ > 0.0) &&
         (binWidthY_ > 0.0)
    )
    {
        label nX = label(rX/binWidthX_);
        label nY = label(rY/binWidthY_);
    
        if
        (
            (nX < nBinsX_) &&
            (nY < nBinsY_) &&
            (nX >= 0) &&
            (nY >= 0)
        )
        {
            inRange = true;
        }
    }

    return inRange;
}



void twoDimForceDistribution::scaleForceDistribution(const scalar& value)
{
    forAll(forces_, x)
    {
        forces_[x] /= value;
        energies_[x] /= value;
    }
}

void twoDimForceDistribution::scaleSampledDistribution()
{
    scaleForceDistribution(mols_);
}

void twoDimForceDistribution::scaleForceDistribution
(
    const List< scalarField>& values
)
{
    //check 
    if(forces_.size() == values.size())
    {
        forAll(forces_, x)
        {
            forAll(forces_[x], y)
            {
                if(values[x][y] > 0.0)
                {
                    forces_[x][y] /= values[x][y];
                    energies_[x][y] /= values[x][y];
                }
            }
        }
    }
    else
    {
        FatalErrorIn("void twoDimForceDistribution::scaleForceDistribution(const scalarField& values)")
            << "Force list is not equal to scaledValues list"
            << exit(FatalError);
    }
}


void twoDimForceDistribution::clear()
{
    //- clear distributions
    forAll(forces_, x)
    {
        forces_[x] = vector::zero;
        mols_[x] = 0.0;
        energies_[x] = 0.0;
    }
}


void twoDimForceDistribution::write
(
    const Time& runTime,
    const polyMesh& mesh
)
{
//     if(runTime.outputTime())
//     {
        //- write out as a lagrangian field

        fileName timePath(runTime.path()/runTime.timeName()/"uniform"/fileName_);

        if (!isDir(timePath))
        {
            mkDir(timePath);
        }

        label nPositions = nBinsX_*nBinsY_;

        //- positions
        OFstream positionsFile(timePath/"positions2D_"+ name_ +".xyz");

        if (positionsFile.good())
        {
            positionsFile << nPositions << endl;
            positionsFile << "(" << endl;

            for (label i = 0; i < nBinsX_; i++)
            {
                for (label j = 0; j < nBinsY_; j++)
                {
                    vector pos = startPoint_ + 
                                 (0.5 + scalar(i))*binWidthX_*unitVectorX_ +
                                 (0.5 + scalar(j))*binWidthY_*unitVectorY_;

                    positionsFile 
                        << "(" << pos.x() << " " << pos.y() << " "
                        << pos.z() << ") " << mesh.findCell(pos)
                        << endl;  
                }
            }

            positionsFile << ")" << endl;
        }
        else
        {
            FatalErrorIn( "twoDimForceDistrPBC::writeField")
                << "Cannot open file "
                << positionsFile.name()
                << abort(FatalError);
        }

        //- positions 2
        OFstream positionsFile2(timePath/"pos2D_"+ name_ +".xyz");

        if (positionsFile2.good())
        {
            positionsFile2 << nPositions << endl;
            positionsFile2 << "(" << endl;

            for (label i = 0; i < nBinsX_; i++)
            {
                for (label j = 0; j < nBinsY_; j++)
                {
                    vector pos = startPoint_ + 

                                (0.5 + scalar(i))*binWidthX_*unitVectorX_ +
                                (0.5 + scalar(j))*binWidthY_*unitVectorY_;

                    positionsFile2 
                        << "(" << pos.x() << " " << pos.y() << " "
                        << pos.z() << ")" 
                        << endl;  
                }
            }

            positionsFile2 << ")" << endl;
        }
        else
        {
            FatalErrorIn( "threeDimForceDistrPBC::writeField")
                << "Cannot open file "
                << positionsFile2.name()
                << abort(FatalError);
        }

//         Info <<"forces: " << nl << forces_ << endl;


        //- forces
        OFstream forcesFile(timePath/"forces2D_"+ name_ +".xyz");

        if (forcesFile.good())
        {
            forcesFile << nPositions << endl;
            forcesFile << "(" << endl;

            for (label i = 0; i < nBinsX_; i++)
            {
                for (label j = 0; j < nBinsY_; j++)
                {
                    forcesFile 
                        << "(" << forces_[i][j].x() << " " << forces_[i][j].y() << " "
                        << forces_[i][j].z() << ") " 
                        << endl;  
                }
            }

            forcesFile << ")" << endl;
        }
        else
        {
            FatalErrorIn( "twoDimForceDistrPBC::writeField")
                << "Cannot open file "
                << forcesFile.name()
                << abort(FatalError);
        }

        //- energies
        OFstream energiesFile(timePath/"energies2D_"+ name_ +".xyz");

        if (energiesFile.good())
        {
            energiesFile << nPositions << endl;
            energiesFile << "(" << endl;

            for (label i = 0; i < nBinsX_; i++)
            {
                for (label j = 0; j < nBinsY_; j++)
                {
                    energiesFile 
                        << energies_[i][j]
                        << endl;
                }
            }

            energiesFile << ")" << endl;
        }
        else
        {
            FatalErrorIn( "twoDimForceDistrPBC::writeField")
                << "Cannot open file "
                << energiesFile.name()
                << abort(FatalError);
        }
//     }
}






// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
