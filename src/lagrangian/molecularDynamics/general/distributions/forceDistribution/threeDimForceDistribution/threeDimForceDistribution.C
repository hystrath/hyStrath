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
    threeDimForceDistribution

Description

\*----------------------------------------------------------------------------*/

#include "threeDimForceDistribution.H"
#include "graph.H"
#include "IFstream.H"
#include "OFstream.H"

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //
// void threeDimForceDistribution::setRadii()
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
threeDimForceDistribution::threeDimForceDistribution()
:
    name_(),
    fileName_("forceDistributions"),
    startPoint_(vector::zero),
    lengthX_(-1.0),
    lengthY_(-1.0),
    lengthZ_(-1.0),
    unitVectorX_(vector::zero),
    unitVectorY_(vector::zero),
    unitVectorZ_(vector::zero),
    nBinsX_(-1),
    nBinsY_(-1),
    nBinsZ_(-1),
    binWidthX_(0.0),
    binWidthY_(0.0),
    binWidthZ_(0.0),
    mols_(),
    forces_(),
    energies_()
{}





// Construct from dictionary
threeDimForceDistribution::threeDimForceDistribution
(
    const dictionary& dict
)
:
    name_(dict.lookup("distributionName")),
    fileName_("forceDistributions"),
    startPoint_(dict.lookup("startPoint")),
    lengthX_(readScalar(dict.lookup("lengthX"))),
    lengthY_(readScalar(dict.lookup("lengthY"))),
    lengthZ_(readScalar(dict.lookup("lengthZ"))),
    unitVectorX_(dict.lookup("unitVectorX")),
    unitVectorY_(dict.lookup("unitVectorY")),
    unitVectorZ_(dict.lookup("unitVectorZ")),
    nBinsX_(readLabel(dict.lookup("nBinsX"))),
    nBinsY_(readLabel(dict.lookup("nBinsY"))),
    nBinsZ_(readLabel(dict.lookup("nBinsZ"))),
    binWidthX_(0.0),
    binWidthY_(0.0),
    binWidthZ_(0.0),
    mols_(nBinsX_),
    forces_(nBinsX_),
    energies_(nBinsX_)
{

    unitVectorX_ /= mag(unitVectorX_);
    unitVectorY_ /= mag(unitVectorY_);
    unitVectorZ_ /= mag(unitVectorZ_);

    binWidthX_ = lengthX_/nBinsX_;
    binWidthY_ = lengthY_/nBinsY_;
    binWidthZ_ = lengthZ_/nBinsZ_;

    forAll(mols_, x)
    {
        mols_[x].setSize(nBinsY_);
        forces_[x].setSize(nBinsY_);
        energies_[x].setSize(nBinsY_);

        forAll(mols_[x], y)
        {
            mols_[x][y].setSize(nBinsZ_, 0.0);
            forces_[x][y].setSize(nBinsZ_, vector::zero);
            energies_[x][y].setSize(nBinsZ_, 0.0);
        }
    }

//     setRadii();
}


// Construct as copy
// threeDimForceDistribution::threeDimForceDistribution(const threeDimForceDistribution& d)
// :
//     Map<label>(static_cast< Map<label> >(d)),
//     binWidth_(d.binWidth())
// {}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

threeDimForceDistribution::~threeDimForceDistribution()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

//- used after null constructor to set properly the data members


void threeDimForceDistribution::addToDistribution
(
    const vector& r, 
    const vector& force,
    const scalar& energy
)
{
    vector rLocal = r - startPoint_;

    scalar rX = rLocal & unitVectorX_;
    scalar rY = rLocal & unitVectorY_;
    scalar rZ = rLocal & unitVectorZ_;

    if
    (
         (binWidthX_ > 0.0) &&
         (binWidthY_ > 0.0) &&
         (binWidthZ_ > 0.0)
    )
    {
        label nX = label(rX/binWidthX_);
        label nY = label(rY/binWidthY_);
        label nZ = label(rZ/binWidthZ_);
    
        if
        (
            (nX < nBinsX_) &&
            (nY < nBinsY_) &&
            (nZ < nBinsZ_) &&
            (nX >= 0) &&
            (nY >= 0) &&
            (nZ >= 0)
        )
        {
            forces_[nX][nY][nZ] += force;
            energies_[nX][nY][nZ] += energy;
            mols_[nX][nY][nZ] += 1.0;
        }
    }
}

bool threeDimForceDistribution::isWithinDistributionRange(const vector& r)
{

    bool inRange = false;

    vector rLocal = r - startPoint_;

    scalar rX = rLocal & unitVectorX_;
    scalar rY = rLocal & unitVectorY_;
    scalar rZ = rLocal & unitVectorZ_;

    if
    (
         (binWidthX_ > 0.0) &&
         (binWidthY_ > 0.0) &&
         (binWidthZ_ > 0.0)
    )
    {
        label nX = label(rX/binWidthX_);
        label nY = label(rY/binWidthY_);
        label nZ = label(rZ/binWidthZ_);
    
        if
        (
            (nX < nBinsX_) &&
            (nY < nBinsY_) &&
            (nZ < nBinsZ_) &&
            (nX >= 0) &&
            (nY >= 0) &&
            (nZ >= 0)
        )
        {
            inRange = true;
        }
    }
    return inRange;
}



void threeDimForceDistribution::scaleForceDistribution(const scalar& value)
{
    forAll(forces_, x)
    {
        forAll(forces_, y)
        {
            forces_[x][y] /= value;
            energies_[x][y] /= value;
        }
    }
}

void threeDimForceDistribution::scaleSampledDistribution()
{
    scaleForceDistribution(mols_);
}

void threeDimForceDistribution::scaleForceDistribution
(
    const List< List< scalarField> >& values
)
{
    //check 
    if(forces_.size() == values.size())
    {
        forAll(forces_, x)
        {
            forAll(forces_[x], y)
            {
                forAll(forces_[x][y], z)
                {
                    if(values[x][y][z] > 0.0)
                    {
                        forces_[x][y][z] /= values[x][y][z];
                        energies_[x][y][z] /= values[x][y][z];
                    }
                }
            }
        }
    }
    else
    {
        FatalErrorIn("void threeDimForceDistribution::scaleForceDistribution(const scalarField& values)")
            << "Force list is not equal to scaledValues list"
            << exit(FatalError);
    }
}

//- return the force distribution in the direction of the vector
//  between the starting and end point

// List< Pair<scalar> > threeDimForceDistribution::fMag()
// {
//     List< Pair<scalar> > forceDistrib(magForces_.size());
// 
//     forAll(forceDistrib, bin)
//     {
//         forceDistrib[bin].first() = magRadii_[bin];
//         forceDistrib[bin].second() = magForces_[bin];
//     }
// 
//     return forceDistrib;
// }


void threeDimForceDistribution::clear()
{
    //- clear distributions
    forAll(forces_, x)
    {
        forAll(forces_[x], y)
        {
            forces_[x][y] = vector::zero;
            mols_[x][y] = 0.0;
            energies_[x][y] = 0.0;
        }
    }
}


void threeDimForceDistribution::write
(
    const Time& runTime,
    const polyMesh& mesh
)
{
//     if(runTime.outputTime())
//     {
        // - rescale the distribution 
//         forAll(forces_, x)
//         {
//             forAll(forces_[x], y)
//             {
//                 forAll(forces_[x][y], z)
//                 {
// 
//                     if(mols_[x][y][z] > 0.0)
//                     {
//                         forces_[x][y][z] /= mols_[x][y][z];
//                         energies_[x][y][z] /= mols_[x][y][z];
//                     }
//                 }       
//             }
//         }

        //- write out as a lagrangian field

        fileName timePath(runTime.path()/runTime.timeName()/"uniform"/fileName_);

        if (!isDir(timePath))
        {
            mkDir(timePath);
        }

        label nPositions = nBinsX_*nBinsY_*nBinsZ_;

        //- positions
        OFstream positionsFile(timePath/"positions3D_"+ name_ +".xyz");

        if (positionsFile.good())
        {
            positionsFile << nPositions << endl;
            positionsFile << "(" << endl;

            for (label i = 0; i < nBinsX_; i++)
            {
                for (label j = 0; j < nBinsY_; j++)
                {
                    for (label k = 0; k < nBinsZ_; k++)
                    {
                        vector pos = startPoint_ + 
                                    (0.5+scalar(i))*binWidthX_*unitVectorX_ +
                                    (0.5+scalar(j))*binWidthY_*unitVectorY_ +
                                    (0.5+scalar(k))*binWidthZ_*unitVectorZ_;
    
                        positionsFile 
                            << "(" << pos.x() << " " << pos.y() << " "
                            << pos.z() << ") " << mesh.findCell(pos)
                            << endl;  
                    }
                }
            }

            positionsFile << ")" << endl;
        }
        else
        {
            FatalErrorIn( "threeDimForceDistrPBC::writeField")
                << "Cannot open file "
                << positionsFile.name()
                << abort(FatalError);
        }

        //- positions 2
        OFstream positionsFile2(timePath/"pos3D_"+ name_ +".xyz");

        if (positionsFile2.good())
        {
            positionsFile2 << nPositions << endl;
            positionsFile2 << "(" << endl;

            for (label i = 0; i < nBinsX_; i++)
            {
                for (label j = 0; j < nBinsY_; j++)
                {
                    for (label k = 0; k < nBinsZ_; k++)
                    {
                        vector pos = startPoint_ + 
                                    (0.5+scalar(i))*binWidthX_*unitVectorX_ +
                                    (0.5+scalar(j))*binWidthY_*unitVectorY_ +
                                    (0.5+scalar(k))*binWidthZ_*unitVectorZ_;
    
                        positionsFile2 
                            << "(" << pos.x() << " " << pos.y() << " "
                            << pos.z() << ")" 
                            << endl;  
                    }
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


        //- forces
        OFstream forcesFile(timePath/"forces3D_"+ name_ +".xyz");

        if (forcesFile.good())
        {
            forcesFile << nPositions << endl;
            forcesFile << "(" << endl;

            for (label i = 0; i < nBinsX_; i++)
            {
                for (label j = 0; j < nBinsY_; j++)
                {
                    for (label k = 0; k < nBinsZ_; k++)
                    {
                        forcesFile 
                            << "(" << forces_[i][j][k].x() << " " << forces_[i][j][k].y() << " "
                            << forces_[i][j][k].z() << ") " 
                            << endl;
                    }
                }
            }

            forcesFile << ")" << endl;
        }
        else
        {
            FatalErrorIn( "threeDimForceDistrPBC::writeField")
                << "Cannot open file "
                << forcesFile.name()
                << abort(FatalError);
        }

        //- energies
        OFstream energiesFile(timePath/"energies3D_"+ name_ +".xyz");

        if (energiesFile.good())
        {
            energiesFile << nPositions << endl;
            energiesFile << "(" << endl;

            for (label i = 0; i < nBinsX_; i++)
            {
                for (label j = 0; j < nBinsY_; j++)
                {
                    for (label k = 0; k < nBinsZ_; k++)
                    {
                        energiesFile 
                            << energies_[i][j][k]  
                            << endl;
                    }
                }
            }

            energiesFile << ")" << endl;
        }
        else
        {
            FatalErrorIn( "threeDimForceDistrPBC::writeField")
                << "Cannot open file "
                << energiesFile.name()
                << abort(FatalError);
        }
//     }
}






// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
