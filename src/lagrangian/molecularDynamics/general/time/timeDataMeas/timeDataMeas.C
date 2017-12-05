/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2007 OpenCFD Ltd.
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

\*---------------------------------------------------------------------------*/

#include "timeDataMeas.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Null Constructor
timeDataMeas::timeDataMeas
(
    Time& t
)
:
    time_(t),
    timeDict_(),
    timeMeasOption_(),
    writeInterval_(readScalar(t.controlDict().lookup("writeInterval"))),
    writeIntSteps_(label((writeInterval_/t.deltaT().value())  + 0.5)),
    resetFieldsAtOutput_(true),
    resetIndex_(0),
    mdTime_(1),
    samplingTime_(),
    averagingTime_(),
    writeTime_(),
    nAvTimeSteps_("nAvTimeSteps_", dimless, 0.0),
    totalNSampSteps_(0),
    totalNAvSteps_(0),
    averagingTimeIndex_(0),
    samplingTimeIndex_(0)
//     averagingTimes_(),
//     samplingTimes_()
{}


//- Construct from Time and timeDict 
timeDataMeas::timeDataMeas
(
    Time& t,
    const dictionary& timeDict
)
:
    time_(t),
    timeDict_(timeDict),
    timeMeasOption_(timeDict_.lookup("timeOption")),
    writeInterval_(readScalar(t.controlDict().lookup("writeInterval"))),
    writeIntSteps_(label((writeInterval_/t.deltaT().value()) + 0.5)),
    resetFieldsAtOutput_(true),
    resetIndex_(0),
    mdTime_(1),
    samplingTime_(),
    averagingTime_(),
    writeTime_(),
    nAvTimeSteps_("nAvTimeSteps_", dimless, 0.0),
    totalNSampSteps_(0),
    totalNAvSteps_(0),
    averagingTimeIndex_(0),
    samplingTimeIndex_(0)
//     averagingTimes_(),
//     samplingTimes_()
{
    setInitialData();
}

void timeDataMeas::setInitialData()
{
    Info << nl << "TimeData Statistics: " << endl;

    scalar deltaTMD = readScalar(time_.controlDict().lookup("deltaT"));

    mdTime_.deltaT() = deltaTMD;

    if (timeDict_.found("resetAtOutput"))
    {
        resetFieldsAtOutput_ = Switch(timeDict_.lookup("resetAtOutput"));
    }

    if(timeMeasOption_ == "write")
    {
        samplingTime_.nSteps() = 1;
        averagingTime_.nSteps() = writeIntSteps_;
    }
    else if(timeMeasOption_ == "decoupledWrite")
    {
        samplingTime_.nSteps() = 1;

        IOdictionary couplingDict
        (
            IOobject
            (
                "couplingDict",
                time().system(),
                time().db(),
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        );

        dictionary timeCouplingDict(couplingDict.subDict("timeCouplingProperties"));

        const scalar couplingTimeInterval = readScalar(timeCouplingDict.lookup("couplingTimeInterval"));
        const scalar molecularOnInterval = readScalar(timeCouplingDict.lookup("molecularOnInterval"));

        scalar ratio = writeInterval_/couplingTimeInterval;
        label nAveragingSteps = label((ratio*molecularOnInterval)/deltaTMD);

        averagingTime_.nSteps() = nAveragingSteps;
    }
    else if(timeMeasOption_ == "general")
    {
        const label nSamples = readLabel(timeDict_.lookup("nSamples"));
        const label nAverages = readLabel(timeDict_.lookup("nAverages"));

        samplingTime_.nSteps() = nSamples;
        averagingTime_.nSteps() = nAverages;

        checkAndModifyTimeProperties();
    }
    else if(timeMeasOption_ == "coupling")
    {
        const label nSamplesDict = readLabel(timeDict_.lookup("nSamples"));

        IOdictionary couplingDict
        (
            IOobject
            (
                "couplingDict",
                time().system(),
                time().db(),
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        );

        dictionary timeCouplingDict(couplingDict.subDict("timeCouplingProperties"));

//         const scalar averagingTime = readScalar(timeCouplingDict.lookup("couplingTimeInterval"));
        const scalar averagingTime = readScalar(timeCouplingDict.lookup("molecularCouplingTimeInterval"));

        samplingTime_.nSteps() = nSamplesDict;
        averagingTime_.nSteps() = label((averagingTime/deltaTMD) + 0.5);

        Info << "Averaging time: " << averagingTime << ", deltaTMD: " << deltaTMD << averagingTime_.nSteps() << endl; 
//         checkAndModifyTimeProperties();

        //- test properties

        bool changedProperties = false;

        const label& nAverages = averagingTime_.nSteps();
        label& nSamples = samplingTime_.nSteps();

        if(nSamples < 1)
        {
            nSamples = 1;
            changedProperties = true;
        }
        else
        {
            if(nSamples > nAverages)
            {
                nSamples = nAverages;
                changedProperties = true;
            }
            else
            {
                while((nAverages % nSamples) != 0)
                {
                    nSamples--;
                    changedProperties = true;
                }
            }
        }

        if(changedProperties)
        {
            Info << "initial nSamples: " << nSamplesDict
                << ", modified nSamples: " << samplingTime_.nSteps() << endl;

            FatalErrorIn("timeDataMeas::timeDataMeas()")
                << "Time properties are inconsistent."
                << " Check and change them appropriately from the time dictionary"
                << nl
                << exit(FatalError);
        }
    }

    else if(timeMeasOption_ == "decoupled")
    {
        const label nSamples = readLabel(timeDict_.lookup("nSamples"));

        IOdictionary couplingDict
        (
            IOobject
            (
                "couplingDict",
                time().system(),
                time().db(),
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        );

        dictionary timeCouplingDict(couplingDict.subDict("timeCouplingProperties"));

        const scalar averagingTime = readScalar(timeCouplingDict.lookup("molecularOnInterval"));

        samplingTime_.nSteps() = nSamples;
        averagingTime_.nSteps() = label(averagingTime/deltaTMD);

        checkAndModifyTimeProperties();
    }
    else if(timeMeasOption_ == "decoupledFromWriteInterval")
    {
        const label nSamplesDict = readLabel(timeDict_.lookup("nSamples"));
        const label nAveragesDict = readLabel(timeDict_.lookup("nAverages"));

        samplingTime_.nSteps() = nSamplesDict;
        averagingTime_.nSteps() = nAveragesDict;

        //- test properties

        bool changedProperties = false;

        const label& nAverages = averagingTime_.nSteps();
        label& nSamples = samplingTime_.nSteps();

        if(nSamples < 1)
        {
            nSamples = 1;
            changedProperties = true;
        }
        else
        {
            if(nSamples > nAverages)
            {
                nSamples = nAverages;
                changedProperties = true;
            }
            else
            {
                while((nAverages % nSamples) != 0)
                {
                    nSamples--;
                    changedProperties = true;
                }
            }
        }

        if(changedProperties)
        {
            Info << "initial nSamples: " << nSamplesDict
                << ", modified nSamples: " << samplingTime_.nSteps() << endl;

            FatalErrorIn("timeDataMeas::timeDataMeas()")
                << "Time properties are inconsistent."
                << " Check and change them appropriately from the time dictionary"
                << nl
                << exit(FatalError);
        }
    }

    else if(timeMeasOption_ == "numberWrite")
    {
        const label noWriteIntervals = readLabel(timeDict_.lookup("noOfWriteIntervals"));
        const label nSamplesDict = readLabel(timeDict_.lookup("nSamples"));

        averagingTime_.nSteps() = (scalar(noWriteIntervals)*writeInterval_)/mdTime_.deltaT();
        samplingTime_.nSteps() = nSamplesDict;
        writeTime_.nSteps() = averagingTime_.nSteps();

        //- test properties

        bool changedProperties = false;

        const label& nAverages = averagingTime_.nSteps();
        label& nSamples = samplingTime_.nSteps();

        if(nSamples < 1)
        {
            nSamples = 1;
            changedProperties = true;
        }
        else
        {
            if(nSamples > nAverages)
            {
                nSamples = nAverages;
                changedProperties = true;
            }
            else
            {
                while((nAverages % nSamples) != 0)
                {
                    nSamples--;
                    changedProperties = true;
                }
            }
        }

        if(changedProperties)
        {
            Info << "initial nSamples: " << nSamplesDict
                << ", modified nSamples: " << samplingTime_.nSteps() << endl;

            FatalErrorIn("timeDataMeas::timeDataMeas()")
                << "Time properties are inconsistent."
                << " Check and change them appropriately from the time dictionary"
                << nl
                << exit(FatalError);
        }
    }
    else
    {
        FatalErrorIn("timeDataMeas::setInitialData()")
            << "timeOption: \"" << timeMeasOption_
            << "\" is not one of the available options." << nl 
            << exit(FatalError);
    }

    samplingTime_.deltaT() = deltaTMD * scalar(samplingTime_.nSteps());
    averagingTime_.deltaT() = deltaTMD * scalar(averagingTime_.nSteps());

    Info << " measurement option: " << timeMeasOption_ << endl;
    Info << " nSamples: " << samplingTime_.nSteps()
         << ", time interval: " << samplingTime_.deltaT()
         << endl;

    Info << " nAverages: " << averagingTime_.nSteps()
         << ", time interval: " << averagingTime_.deltaT()
         << endl;


    const scalar& endTime = time_.endTime().value();
    const scalar& startTime = time_.startTime().value();

    totalNAvSteps_ = label ( ((endTime - startTime) / averagingTime_.deltaT()) + 0.5 );

//     Info<< "test: " << " av time: " << averagingTime_.deltaT() 
//         << " totalNavSteps: " << (endTime - startTime) / averagingTime_.deltaT() 
//         << " label: " << label((endTime - startTime) / averagingTime_.deltaT())
//         << " adjustement : " 
//         << label( ( (endTime - startTime) / averagingTime_.deltaT() ) + 0.5 )
//         << " final result : " << totalNAvSteps_
//         << endl;

//     totalNCalcSteps_ = label((endTime - startTime) / calcPropTime_.deltaT());

    totalNSampSteps_ = label(((endTime - startTime) / samplingTime_.deltaT()) + 0.5);

    Info << " total no. of sampling steps: " << totalNSampSteps_ << endl;
    Info << " total no. of averaging Steps: " << totalNAvSteps_ << endl;


    Info << nl << endl;

//     averagingTimes_.setSize(totalNAvSteps_+1, 0.0);
// 
//     forAll(averagingTimes_, tT)
//     {
//         averagingTimes_[tT] = startTime + tT*averagingTime_.deltaT();
//     }

    nAvTimeSteps_.value() = scalar(averagingTime_.nSteps());

//     samplingTimes_.setSize(totalNSampSteps_+1, 0.0);
// 
//     forAll(samplingTimes_, tT)
//     {
//         samplingTimes_[tT] = startTime + tT*samplingTime_.deltaT();
//     }

}


void timeDataMeas::checkAndModifyTimeProperties()
{
    //- checking 

    bool changedProperties = false;

    const label nAveragesOriginal = averagingTime_.nSteps();
    const label nSamplesOriginal = samplingTime_.nSteps();

    // - 1. averaging time

    label& nAverages = averagingTime_.nSteps();

    if(nAverages < 1)
    {
        nAverages = 1;
        changedProperties = true;
    }
    else
    {
        if(nAverages > writeIntSteps_)
        {
            nAverages = writeIntSteps_;
            changedProperties = true;
        }
        else
        {
            while((writeIntSteps_ % nAverages) != 0)
            {
                nAverages++;
                changedProperties = true;
            }
        }
    }

    // - 2. sampling time
    label& nSamples = samplingTime_.nSteps();

    if(nSamples < 1)
    {
        nSamples = 1;
        changedProperties = true;
    }
    else
    {
        if(nSamples > nAverages)
        {
            nSamples = nAverages;
            changedProperties = true;
        }
        else
        {
            while((nAverages % nSamples) != 0)
            {
                nSamples--;
                changedProperties = true;
            }
        }
    }

    if(changedProperties)
    {
        Info << "initial nSamples: " << nSamplesOriginal
             << ", modified nSamples: " << samplingTime_.nSteps() << endl;

        Info << "initial nAverages: " << nAveragesOriginal
             << ", modified nAverages: " << averagingTime_.nSteps() << endl;

        FatalErrorIn("timeDataMeas::timeDataMeas()")
            << "Time properties are inconsistent."
            << " Check and change them appropriately from the time dictionary" 
            << nl 
            << exit(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

timeDataMeas::~timeDataMeas()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void timeDataMeas::setTimeData(const dictionary& timeDict)
{
    const label nSamples(readLabel(timeDict.lookup("nSamples")));
    const label nAverages(readLabel(timeDict.lookup("nAverages")));
//     const label nCalcProp(readLabel(timeDict.lookup("nCalcProp")));
    
 
    samplingTime_.nSteps() = nSamples;
    averagingTime_.nSteps() = nAverages;
//     calcPropTime_.nSteps() = nCalcProp;
    

    setInitialData();
}


const Time& timeDataMeas::time() const
{
    return time_;
}

Time& timeDataMeas::time()
{
    return time_;
}


const scalar& timeDataMeas::writeInterval() const
{
    return writeInterval_;
}

const label& timeDataMeas::writeIntervalSteps() const
{
    return writeIntSteps_;
}

const bool& timeDataMeas::samplingTime() const
{
    return samplingTime_.endTime();
}

const bool& timeDataMeas::averagingTime() const
{
    return averagingTime_.endTime();
}

const bool& timeDataMeas::writeTime() const
{
    return writeTime_.endTime();
}

// const bool& timeDataMeas::calcPropTime() const
// {
//     return calcPropTime_.endTime();
// }
// 


const label& timeDataMeas::nSamples() const
{
    return samplingTime_.nSteps();
}


const label& timeDataMeas::nAverages() const
{
    return averagingTime_.nSteps();
}

const dimensionedScalar& timeDataMeas::nAvTimeSteps() const
{
    return nAvTimeSteps_;
}

// for resetting
scalar timeDataMeas::nAveragingTimeSteps()
{
    return scalar(nAvTimeSteps().value()*resetIndex_);
}
// const label& timeDataMeas::nCalcProp() const
// {
//     return calcPropTime_.nSteps();
// }

const label& timeDataMeas::totalNSampSteps() const
{
    return totalNSampSteps_;
}

const label& timeDataMeas::totalNAvSteps() const
{
    return totalNAvSteps_;
}

const label& timeDataMeas::averagingTimeIndex() const
{
    return averagingTimeIndex_;
}

const label& timeDataMeas::samplingTimeIndex() const
{
    return samplingTimeIndex_;
}


const bool& timeDataMeas::resetFieldsAtOutput() const
{
    return resetFieldsAtOutput_;
}

bool& timeDataMeas::resetFieldsAtOutput()
{
    return resetFieldsAtOutput_;
}




// const label& timeDataMeas::totalNCalcSteps() const
// {
//     return totalNCalcSteps_;
// }

// const label& timeDataMeas::calcTimeIndex() const
// {
//     return calcTimeIndex_;
// }


/*
const scalarField& timeDataMeas::averagingTimes() const
{
    return averagingTimes_;
}
const scalarField& timeDataMeas::samplingTimes() const
{
    return samplingTimes_;
}
*/

scalarField timeDataMeas::averagingTimes()
{
//     const scalar& endTime = time_.endTime().value();
    const scalar& startTime = time_.startTime().value();

    scalarField averagingTimes(totalNAvSteps_+1, 0.0);

    forAll(averagingTimes, tT)
    {
        averagingTimes[tT] = startTime + tT*averagingTime_.deltaT();
    }

    return averagingTimes;
}

scalarField timeDataMeas::samplingTimes()
{
//     const scalar& endTime = time_.endTime().value();
    const scalar& startTime = time_.startTime().value();

    scalarField samplingTimes(totalNSampSteps_+1, 0.0);

    forAll(samplingTimes, tT)
    {
        samplingTimes[tT] = startTime + tT*samplingTime_.deltaT();
    }

    return samplingTimes;
}

scalarField timeDataMeas::writeTimes()
{
    const scalar& endTime = time_.endTime().value();
    const scalar& startTime = time_.startTime().value();

    label nWriteSteps = label((endTime - startTime)/writeInterval_);

    scalarField writingTimes(nWriteSteps+1, 0.0);

    forAll(writingTimes, tT)
    {
        writingTimes[tT] = startTime + tT*writeInterval_;
    }

    return writingTimes;
}


label timeDataMeas::writeSteps()
{
    const scalar& endTime = time_.endTime().value();
    const scalar& startTime = time_.startTime().value();

    label nWriteSteps = label((endTime - startTime)/writeInterval_);

    return nWriteSteps;
}



const timeInterval& timeDataMeas::mdTimeInterval() const
{
    return mdTime_;
}

const timeInterval& timeDataMeas::sampleTimeInterval() const
{
    return samplingTime_;
}

const timeInterval& timeDataMeas::averageTimeInterval() const
{
    return averagingTime_;
}

const timeInterval& timeDataMeas::writeTimeInterval() const
{
    return writeTime_;
}


const word& timeDataMeas::timeOption() const
{
    return timeMeasOption_;
}



// const timeInterval& timeDataMeas::calcPropTimeInterval() const
// {
//     return calcPropTime_;
// }

// timeInterval& timeDataMeas::calcPropTimeInterval()
// {
//     return calcPropTime_;
// }


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

//- Prefix increment
timeDataMeas& timeDataMeas::operator++()
{
    samplingTime_++;
    averagingTime_++;
    writeTime_++;

//     calcPropTime_++;
    if(averagingTime_.endTime())
    {

        if(resetFieldsAtOutput_)
        {
            resetIndex_ = -1;
        }

        averagingTimeIndex_++;
        resetIndex_++;
    }
    if(samplingTime_.endTime())
    {
        samplingTimeIndex_++;
    }
//     if(calcPropTime_.endTime())
//     {
//         calcTimeIndex_++;
//     }


    return *this;
}


//- Postfix increment
timeDataMeas& timeDataMeas::operator++(int)
{
    return operator++();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
