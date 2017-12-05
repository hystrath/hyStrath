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
    forceDistribution

Description

\*----------------------------------------------------------------------------*/

#include "forceDistribution.H"
#include "graph.H"
#include "IFstream.H"
#include "OFstream.H"

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //
void forceDistribution::setRadii()
{
    for(label i = 0; i < noOfBins_; i++)
    {
        magRadii_[i] = 0.5*binWidth_ + scalar(i)*binWidth_;

        radii_[i] = startPoint_ + (0.5 + scalar(i))*binWidth_*unitVector_;
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Null constructor
forceDistribution::forceDistribution()
:
    name_(),
    fileName_("forceDistributions"),
    startPoint_(vector::zero),
    endPoint_(vector::zero),
    unitVector_(vector::zero),
    noOfBins_(-1),
    forces_(),
    energies_(),
    virials_(),
    mols_(),
    radii_(),
    magRadii_(),
    magForces_(),
    binWidth_(0.0)
{}

// Construct from name --- used for reading a distribution from file
forceDistribution::forceDistribution
(
    const word& name
//     const word& fileName
)
:
    name_(name),
    fileName_("forceDistributions"),
    startPoint_(vector::zero),
    endPoint_(vector::zero),
    unitVector_(vector::zero),
    noOfBins_(-1),
    forces_(),
    energies_(),
    virials_(),
    mols_(),
    radii_(),
    magRadii_(),
    magForces_(),
    binWidth_(0.0)
{}


// Construct from components --- initialises properly the data members
forceDistribution::forceDistribution
(
    const word& name,
//     const word& fileName,
    const vector& startPoint,
    const vector& endPoint,
    const label& noOfBins
)
:
    name_(name),
    fileName_("forceDistributions"),
    startPoint_(startPoint),
    endPoint_(endPoint),
    unitVector_
    (
        (endPoint_- startPoint_)/
        mag(endPoint_- startPoint_)
    ),
    noOfBins_(noOfBins),
    forces_(noOfBins, vector::zero),
    energies_(noOfBins, scalar(0.0)),
    virials_(noOfBins, scalar(0.0)),
    mols_(noOfBins, scalar(0.0)),
    radii_(noOfBins, vector::zero),
    magRadii_(noOfBins, scalar(0.0)),
    magForces_(noOfBins, scalar(0.0)),
    binWidth_(mag(endPoint_ - startPoint_)/(noOfBins_))
{
    setRadii();
}

// Construct from dictionary
forceDistribution::forceDistribution
(
    const dictionary& dict
)
:
    name_(dict.lookup("distributionName")),
//     fileName_(dict.lookup("fileName")),
    fileName_("forceDistributions"),
    startPoint_(dict.lookup("startPoint")),
    endPoint_(dict.lookup("endPoint")),
    unitVector_
    (
        (endPoint_- startPoint_)/
        mag(endPoint_- startPoint_)
    ),
    noOfBins_(readLabel(dict.lookup("noOfBins"))),
    forces_(noOfBins_, vector::zero),
    energies_(noOfBins_, scalar(0.0)),
    virials_(noOfBins_, scalar(0.0)),
    mols_(noOfBins_, scalar(0.0)),
    radii_(noOfBins_, vector::zero),
    magRadii_(noOfBins_, scalar(0.0)),
    magForces_(noOfBins_, scalar(0.0)),
    binWidth_(mag(endPoint_- startPoint_)/(noOfBins_))
{
    setRadii();
}


// Construct as copy
// forceDistribution::forceDistribution(const forceDistribution& d)
// :
//     Map<label>(static_cast< Map<label> >(d)),
//     binWidth_(d.binWidth())
// {}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

forceDistribution::~forceDistribution()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

//- used after null constructor to set properly the data members
void forceDistribution::setForceDistr
(
    const word& name,
//     const word& fileName,
    const vector& startPoint,
    const vector& endPoint,
    const label& noOfBins
)
{
    name_ = name;
//     fileName_ = fileName;
    startPoint_ = startPoint;
    endPoint_ = endPoint;
    unitVector_ = (endPoint_- startPoint_)/mag(endPoint_- startPoint_);
    noOfBins_ = noOfBins;
    binWidth_ = mag(endPoint_ - startPoint_)/(noOfBins_);

    forces_.setSize(noOfBins, vector::zero);
    energies_.setSize(noOfBins, scalar(0.0));
    virials_.setSize(noOfBins, scalar(0.0));
    mols_.setSize(noOfBins, scalar(0.0));
    radii_.setSize(noOfBins, vector::zero);
    magRadii_.setSize(noOfBins, scalar(0.0));
    magForces_.setSize(noOfBins, scalar(0.0));

    setRadii();
}

void forceDistribution::setForceDistr
(
    const dictionary& dict
)
{
    const word name(dict.lookup("distributionName"));
    const vector startPoint(dict.lookup("startPoint"));
    const vector endPoint(dict.lookup("endPoint"));
    const label noOfBins(readLabel(dict.lookup("noOfBins")));

    setForceDistr
    (
        name,
        startPoint,
        endPoint,
        noOfBins
    );
}

void forceDistribution::setMagForceDistr()
{
    forAll(magForces_, f)
    {
        if(mag(forces_[f]) > 0.0)
        {
            magForces_[f] = forces_[f] & unitVector_;
        }
    }
}

void forceDistribution::addToDistribution
(
    const vector& r, 
    const vector& force,
    const scalar& energy
)
{
    vector rLocal = r - startPoint_;

    scalar rNormal = rLocal & unitVector_;

    if(binWidth_ > 0.0)
    {
        label n = label(rNormal/binWidth_);
    
        if((n < magRadii_.size()) && (n >= 0))
        {
            forces_[n] += force;
            energies_[n] += energy;
            mols_[n] += 1.0;
        }
    }
}

void forceDistribution::addToDistribution
(
    const vector& r, 
    const vector& force,
    const scalar& energy,
    const scalar& virial
)
{
    vector rLocal = r - startPoint_;

    scalar rNormal = rLocal & unitVector_;

    if(binWidth_ > 0.0)
    {
        label n = label(rNormal/binWidth_);
    
        if((n < magRadii_.size()) && (n >= 0))
        {
            forces_[n] += force;
            energies_[n] += energy;
            virials_[n] += virial;
            mols_[n] += 1.0;
        }
    }
}

bool forceDistribution::isWithinDistributionRange(const vector& r)
{
    vector rLocal = r - startPoint_;

    scalar rNormal = rLocal & unitVector_;

    bool inRange = false;

    if(binWidth_ > 0.0)
    {
        label n = label(rNormal/binWidth_);
    
        if((n < magRadii_.size()) && (n >= 0))
        {
            inRange = true;
        }
    }

    return inRange;
}



void forceDistribution::scaleForceDistribution(const scalar& value)
{
    forces_ /= value;
    energies_ /= value;
    virials_ /= value;
}

void forceDistribution::scaleSampledDistribution()
{
    scaleForceDistribution(mols_);
}

void forceDistribution::scaleForceDistribution
(
    const scalarField& values
)
{
    //check 
    if(forces_.size() == values.size())
    {
        forAll(forces_, f)
        {
            if(values[f] > 0.0)
            {
                forces_[f] /= values[f];

                energies_[f] /= values[f];

                virials_[f] /= values[f];
            }
        }
    }
    else
    {
        FatalErrorIn("void forceDistribution::scaleForceDistribution(const scalarField& values)")
            << "Force list is not equal to scaledValues list"
            << exit(FatalError);
    }
}

//- return the force distribution in the direction of the vector
//  between the starting and end point

List< Pair<scalar> > forceDistribution::fMagDistr()
{
    List< Pair<scalar> > forceDistrib(magForces_.size());

    forAll(forceDistrib, bin)
    {
        forceDistrib[bin].first() = magRadii_[bin];
        forceDistrib[bin].second() = magForces_[bin];
    }

    return forceDistrib;
}

List< Pair<scalar> > forceDistribution::pEpairDistr()
{
    List< Pair<scalar> > pEdistr(energies_.size());

    forAll(pEdistr, bin)
    {
        pEdistr[bin].first() = magRadii_[bin];
        pEdistr[bin].second() = energies_[bin];
    }

    return pEdistr;
}

List< Pair<scalar> > forceDistribution::virialPairDistr()
{
    List< Pair<scalar> > virialDistr(virials_.size());

    forAll(virialDistr, bin)
    {
        virialDistr[bin].first() = magRadii_[bin];
        virialDistr[bin].second() = virials_[bin];
    }

    return virialDistr;
}

void forceDistribution::clear()
{
    forces_ = vector::zero;
    energies_ = 0.0;
    virials_ = 0.0;
    magForces_ = 0.0;
    mols_ = 0.0;
}


void forceDistribution::write(const Time& runTime)
{
//     if(runTime.outputTime())
//     {
    fileName timePath(runTime.path()/runTime.timeName()/"uniform"/fileName_);

    if (!isDir(timePath))
    {
        mkDir(timePath);
    }

    //- write raw file

    setMagForceDistr();

    OFstream magForceDistrFile(timePath/name_+"_force.raw");

    if (magForceDistrFile.good())
    {
        magForceDistrFile << fMagDistr() << endl;
    }
    else
    {
        FatalErrorIn("void forceDistribution::write()")
            << "Cannot open file " << magForceDistrFile.name()
            << abort(FatalError);
    }


    OFstream pEdistrFile(timePath/name_+"_energy.raw");

    if (pEdistrFile.good())
    {
        pEdistrFile << pEpairDistr() << endl;
    }
    else
    {
        FatalErrorIn("void forceDistribution::write()")
            << "Cannot open file " << pEdistrFile.name()
            << abort(FatalError);
    }

    OFstream virialDistrFile(timePath/name_+"_virial.raw");

    if (virialDistrFile.good())
    {
        virialDistrFile << virialPairDistr() << endl;
    }
    else
    {
        FatalErrorIn("void forceDistribution::write()")
            << "Cannot open file " << virialDistrFile.name()
            << abort(FatalError);
    }


    writeTimeData(timePath, name_+"_vectorForce", magRadii_, forces_);
    writeTimeData(timePath, name_+"_magForce", magRadii_, magForces_);
    writeTimeData(timePath, name_+"_energy", magRadii_, energies_);
    writeTimeData(timePath, name_+"_virial", magRadii_, virials_);
//     }
}

//- used just to write the linear distribution which has been read in.
void forceDistribution::writeLinearForceDistribution(const Time& runTime)
{
    fileName timePath(runTime.path()/runTime.timeName()/"uniform"/fileName_);

    if (!isDir(timePath))
    {
        mkDir(timePath);
    }

    OFstream magForceDistrFile(timePath/name_+"_force.raw");

    if (magForceDistrFile.good())
    {
        magForceDistrFile << fMagDistr() << endl;
    }
    else
    {
        FatalErrorIn("void forceDistribution::writeLinearForceDistribution()")
            << "Cannot open file " << magForceDistrFile.name()
            << abort(FatalError);
    }

    OFstream pEdistrFile(timePath/name_+"_energy.raw");

    if (pEdistrFile.good())
    {
        pEdistrFile << pEpairDistr() << endl;
    }
    else
    {
        FatalErrorIn("void forceDistribution::writeLinearForceDistribution()")
            << "Cannot open file " << pEdistrFile.name()
            << abort(FatalError);
    }

    OFstream virialDistrFile(timePath/name_+"_virial.raw");

    if (virialDistrFile.good())
    {
        virialDistrFile << virialPairDistr() << endl;
    }
    else
    {
        FatalErrorIn("void forceDistribution::writeLinearForceDistribution()")
            << "Cannot open file " << virialDistrFile.name()
            << abort(FatalError);
    }

    writeTimeData(timePath, name_+"_magForce", magRadii_, magForces_);
    writeTimeData(timePath, name_+"_energy", magRadii_, energies_);
    writeTimeData(timePath, name_+"_virial", magRadii_, virials_);
}


void forceDistribution::writeTimeData
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


void forceDistribution::writeTimeData
(
    const fileName& pathName,
    const word& nameFile,
    const scalarField& xData,
    const vectorField& yData
)
{
    OFstream file(pathName/nameFile + ".xyz");

    if(file.good())
    {
        forAll(yData, n)
        {
            file 
                << xData[n] << "\t" 
                << yData[n].x() << "\t" << yData[n].y() 
                << "\t" << yData[n].z() 
                << endl;
        }
    }
    else
    {
        FatalErrorIn("void stateController::writeTimeData()")
            << "Cannot open file " << file.name()
            << abort(FatalError);
    }
}


void forceDistribution::read
(
    const Time& runTime
)
{
    fileName timePath(runTime.path()/runTime.timeName()/"uniform"/fileName_);

    if (!isDir(timePath))
    {
        mkDir(timePath);

        FatalErrorIn("void forceDistribution::read()")
            << "Cannot open file " << timePath 
            << abort(FatalError);
    }

    //- read force distribution

    IFstream forceDistrFile(timePath/name_+"_force.raw");

    List< Pair<scalar> > forces;

    if (forceDistrFile.good())
    {
        forceDistrFile >> forces;
    }
    else
    {
        FatalErrorIn("void forceDistribution::read()")
            << "Cannot open file " << forceDistrFile.name()
            << abort(FatalError);
    }

    setForceDistribution(forces);

    // set the bin width of force distribution
    binWidth_ = readBinWidth();
//     Info << "bin Width : " << binWidth_ << endl;


    //-read potential energy
    IFstream energyDistrFile(timePath/name_+"_energy.raw");

    List< Pair<scalar> > energies;

    if (energyDistrFile.good())
    {
        energyDistrFile >> energies;
    }
    else
    {
        FatalErrorIn("void forceDistribution::read()")
            << "Cannot open file " << energyDistrFile.name()
            << abort(FatalError);
    }

    setPEDistribution(energies);

    //-read virial 

    IFstream virialDistrFile(timePath/name_+"_virial.raw");

    List< Pair<scalar> > virials;

    if (virialDistrFile.good())
    {
        virialDistrFile >> virials;
    }
    else
    {
        FatalErrorIn("void forceDistribution::read()")
            << "Cannot open file " << virialDistrFile.name()
            << abort(FatalError);
    }

    setVirialDistribution(virials);

}

void forceDistribution::setForceDistribution
(
    const List< Pair<scalar> >& forces
)
{
    magForces_.setSize(forces.size());
    magRadii_.setSize(forces.size());
    energies_.setSize(forces.size());
    virials_.setSize(forces.size());

    forAll(forces, bin)
    {
        magRadii_[bin] = forces[bin].first();
        magForces_[bin] = forces[bin].second();
    }
}


void forceDistribution::setPEDistribution
(
    const List< Pair<scalar> >& energies
)
{
    if(energies.size() != energies_.size())
    {
        FatalErrorIn("void forceDistribution::setPEDistribution()")
            << "Energy distribution size : "<< energies.size() 
            << " not equal to the force distribution list: " 
            << energies_.size()
            << abort(FatalError);
    }

    forAll(energies, bin)
    {
        energies_[bin] = energies[bin].second();
    }
}

void forceDistribution::setVirialDistribution
(
    const List< Pair<scalar> >& virials
)
{
    if(virials.size() != virials_.size())
    {
        FatalErrorIn("void forceDistribution::setPEDistribution()")
            << "Virial distribution size : "<< virials.size() 
            << " not equal to the force distribution list: " 
            << virials_.size()
            << abort(FatalError);
    }

    forAll(virials, bin)
    {
        virials_[bin] = virials[bin].second();
    }
}



// return bin value for force
scalar forceDistribution::binValue(const scalar& rD) const
{
    label key = label(rD/binWidth_);

    if
    (
        (key < magForces_.size()) &&
        (key >= 0)
    )
    {
        return magForces_[key];
    }
    else 
    {
        return 0.0;
    }
}


scalar forceDistribution::binValuePE(const scalar& rD) const
{
    label key = label(rD/binWidth_);

    if
    (
        (key < energies_.size()) &&
        (key >= 0)
    )
    {
        return energies_[key];
    }
    else 
    {
        return 0.0;
    }
}

scalar forceDistribution::binValueVirial(const scalar& rD) const
{
    label key = label(rD/binWidth_);

    if
    (
        (key < virials_.size()) &&
        (key >= 0)
    )
    {
        return virials_[key];
    }
    else 
    {
        return 0.0;
    }
}

// return linear interpolated value for the force
scalar forceDistribution::linearInterp(const scalar& r) const
{
    label key = label(r/binWidth_);

    if
    (
        (key < magForces_.size()) &&
        (key >= 0)
    )
    {
        scalar deltaR = magRadii_[key] - r;

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

        const scalar& f1 = magForces_[key];
        const scalar& f2 = magForces_[key2];

        scalar fNew = (((f2 - f1)/binWidth_)*mag(deltaR)) + f1;

        return fNew;
    }
    else
    {
        return 0.0;
    }
}


scalar forceDistribution::linearInterpPE(const scalar& r) const
{
    label key = label(r/binWidth_);

    if
    (
        (key < energies_.size()) &&
        (key >= 0)
    )
    {
        scalar deltaR = magRadii_[key] - r;

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

        const scalar& f1 = energies_[key];
        const scalar& f2 = energies_[key2];

        scalar fNew = (((f2 - f1)/binWidth_)*mag(deltaR)) + f1;

        return fNew;
    }
    else
    {
        return 0.0;
    }
}

scalar forceDistribution::linearInterpVirial(const scalar& r) const
{
    label key = label(r/binWidth_);

    if
    (
        (key < virials_.size()) &&
        (key >= 0)
    )
    {
        scalar deltaR = magRadii_[key] - r;

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

        const scalar& f1 = virials_[key];
        const scalar& f2 = virials_[key2];

        scalar fNew = (((f2 - f1)/binWidth_)*mag(deltaR)) + f1;

        return fNew;
    }
    else
    {
        return 0.0;
    }
}

scalar forceDistribution::readBinWidth()
{
    //- test for zero size
    if(magForces_.size() == 0)
    {
        FatalErrorIn("forceDistribution")
            << "Size of forceDistribution:  " << name_
            << " is zero" 
            << abort(FatalError);
    }

    //- test for common bin width

    scalar binWidth = 0.0;

    bool constantBinWidth = true;

    scalar tolerance = 1e-4;

    for(label i = 0; i < magRadii_.size()-1; i++)
    {
        if(binWidth == 0.0)
        {
            binWidth = mag(magRadii_[i+1] - magRadii_[i]);
        }
        else
        {
            scalar newBinWid = mag(magRadii_[i+1] - magRadii_[i]);

            if(mag(newBinWid - binWidth) > tolerance)
            {
                constantBinWidth = false;
                break;
            }
        }
    }

    if(!(binWidth > 0.0) || !constantBinWidth)
    {
        FatalErrorIn("forceDistribution")
            << "Check binWidth in forceDistribution: " << name_
            << " for constant binWidth: " << binWidth
            << abort(FatalError);
    }

    return binWidth;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

// void forceDistribution::operator=(const forceDistribution& rhs)
// {
//     // Check for assignment to self
//     if (this == &rhs)
//     {
//         FatalErrorIn("forceDistribution::operator=(const forceDistribution&)")
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

// Ostream& operator<<(Ostream& os, const forceDistribution& d)
// {
//     os  << d.binWidth_
//         << static_cast<const Map<label>&>(d);
// 
//     // Check state of Ostream
//     os.check
//     (
//         "Ostream& operator<<(Ostream&, "
//         "const forceDistribution&)"
//     );
// 
//     return os;
// }

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
