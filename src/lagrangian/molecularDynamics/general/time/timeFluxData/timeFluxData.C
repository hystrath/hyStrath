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

#include "timeFluxData.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Null Constructor
timeFluxData::timeFluxData
(
    Time& t
)
:
    time_(t),
    writeInterval_(readScalar(t.controlDict().lookup("writeInterval"))),
    writeIntSteps_(label(writeInterval_/t.deltaT().value() + 0.5)),
    mdTime_(1),
    averagingTime_(),
    controlTime_(),
    writeTime_(),
    nAvTimeSteps_("nAvTimeSteps_", dimless, 0.0),
    nControlSteps_(0.0),
    totalNAvSteps_(0),
    totalNContSteps_(0),
    totalNWrSteps_(0),
    controlTimeIndex_(0),
    averagingTimeIndex_(0)
/*    controlTimes_(),
    averagingTimes_()*/
{}


//- Construct from Time and timeDict 
timeFluxData::timeFluxData
(
    Time& t,
    const dictionary& timeDict
)
:
    time_(t),
    writeInterval_(readScalar(t.controlDict().lookup("writeInterval"))),
    writeIntSteps_(label(writeInterval_/t.deltaT().value() + 0.5)),
    mdTime_(1, readScalar(time_.controlDict().lookup("deltaT"))),
    averagingTime_(readLabel(timeDict.lookup("nAverages"))),
    controlTime_(readLabel(timeDict.lookup("nControls"))),
    writeTime_(writeIntSteps_, writeInterval_),
    nAvTimeSteps_("nAvTimeSteps_", dimless, 0.0),
    nControlSteps_(0.0),
    totalNAvSteps_(0),
    totalNContSteps_(0),
    totalNWrSteps_(0),
    controlTimeIndex_(0),
    averagingTimeIndex_(0)
//     controlTimes_(),
//     averagingTimes_()
{
    setInitialData();
}

void timeFluxData::setInitialData()
{
    Info << nl << "TimeData Statistics: " << endl;

//     scalar deltaTMD = readScalar(time_.controlDict().lookup("deltaT"));

//     mdTime_.deltaT() = deltaTMD;

    checkAndModifyTimeProperties();

    const scalar& endTime = time_.endTime().value();
    const scalar& startTime = time_.startTime().value();

    totalNAvSteps_ = label((endTime - startTime) / averagingTime_.deltaT());
    totalNContSteps_ = label((endTime - startTime) / controlTime_.deltaT());
    totalNWrSteps_ = label((endTime - startTime) / writeTime_.deltaT());

    Info << " total no. of averaging Steps: " << totalNAvSteps_ << endl;
    Info << " total no. of control Steps: " << totalNContSteps_ << endl;

    Info << nl << endl;

//     controlTimes_.setSize(totalNContSteps_+1, 0.0);
// 
//     averagingTimes_.setSize(totalNAvSteps_+1, 0.0);
// 
// 
//     forAll(controlTimes_, tT)
//     {
//         controlTimes_[tT] = startTime + tT*controlTime_.deltaT();
//     }
// 
//     forAll(averagingTimes_, tT)
//     {
//         averagingTimes_[tT] = startTime + tT*averagingTime_.deltaT();
//     }

    nAvTimeSteps_.value() = scalar(averagingTime_.nSteps());

    nControlSteps_ = averagingTime_.deltaT()/controlTime_.deltaT();


   //-- offsetting the controlling time index so that the time-interval finishes 
    //   one time-step ahead of the calcProp time.

    averagingTime_.timeIndex() = averagingTime_.nSteps();
    controlTime_.timeIndex()--;
}


void timeFluxData::checkAndModifyTimeProperties()
{
    //- checking 

    bool changedProperties = false;
    const scalar& deltaTMD = mdTime_.deltaT();

    // - 1. averaging time
    //   for now we ensure that the averaging interval is equal 
    //   to the writing interval

    label& nAverages = averagingTime_.nSteps();

    Info << " nAveraging steps (initial): " << nAverages;

    if(nAverages < 1)
    {
        nAverages = 1;
        changedProperties = true;
    }
//     else
//     {
//         if(nAverages > writeIntSteps_)
//         {
//             nAverages = writeIntSteps_;
//             changedProperties = true;
//         }
//         else
//         {
//             while((writeIntSteps_ % nAverages) != 0)
//             {
//                 nAverages++;
//                 changedProperties = true;
//             }
//         }
//     }

    averagingTime_.deltaT() = deltaTMD * scalar(nAverages);

    Info << ", (final): " << nAverages
         << " time interval: " << averagingTime_.deltaT()
         << endl;



   // - 2. control time

    label& nControls = controlTime_.nSteps();

    Info << " nControls (initial): " << nControls;

    if(nControls < 1)
    {
        nControls = 1;
        changedProperties = true;
    }
    else
    {
//         if(nControls > writeIntSteps_)
//         {
//             nControls = writeIntSteps_;
//             changedProperties = true;
//         }

        if(nControls > nAverages)
        {
            nControls = nAverages;
            changedProperties = true;
        }
        else
        {
            while((nAverages % nControls) != 0)
            {
                nControls ++;
                changedProperties = true;
            }
        }
    }

    controlTime_.deltaT() = deltaTMD * scalar(nControls);

    Info << ", (final): " << nControls
         << " time interval: " << controlTime_.deltaT()
         << endl;

    if(changedProperties)
    {
        FatalErrorIn("timeFluxData::timeFluxData()")
            << "Time data members have been changed."
            << " Check and change them appropriately from the time dictionary" 
            << nl 
            << exit(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

timeFluxData::~timeFluxData()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void timeFluxData::setTimeData(const dictionary& timeDict)
{
    const label nAverages(readLabel(timeDict.lookup("nAverages")));
    const label nControls(readLabel(timeDict.lookup("nControls")));
 
    averagingTime_.nSteps() = nAverages;
    controlTime_.nSteps() = nControls;

    setInitialData();
}


const Time& timeFluxData::time() const
{
    return time_;
}

Time& timeFluxData::time()
{
    return time_;
}

const scalar& timeFluxData::writeInterval() const
{
    return writeInterval_;
}

const label& timeFluxData::writeIntervalSteps() const
{
    return writeIntSteps_;
}

const bool& timeFluxData::averagingTime() const
{
    return averagingTime_.endTime();
}

const bool& timeFluxData::controlTime() const
{
    return controlTime_.endTime();
}


const label& timeFluxData::nControls() const
{
    return controlTime_.nSteps();
}

const label& timeFluxData::nAverages() const
{
    return averagingTime_.nSteps();
}

const dimensionedScalar& timeFluxData::nAvTimeSteps() const
{
    return nAvTimeSteps_;
}

const scalar& timeFluxData::nControlSteps() const
{
    return nControlSteps_;
}



const label& timeFluxData::totalNAvSteps() const
{
    return totalNAvSteps_;
}

const label& timeFluxData::totalNContSteps() const
{
    return totalNContSteps_;
}

const label& timeFluxData::totalNWrSteps() const
{
    return totalNWrSteps_;
}




const label& timeFluxData::controlTimeIndex() const
{
    return controlTimeIndex_;
}

const label& timeFluxData::averagingTimeIndex() const
{
    return averagingTimeIndex_;
}

const label& timeFluxData::writeTimeIndex() const
{
    return writeTimeIndex_;
}

// const scalarField& timeFluxData::controlTimes() const
// {
//     return controlTimes_;
// }
// 
// const scalarField& timeFluxData::averagingTimes() const
// {
//     return averagingTimes_;
// }


scalarField timeFluxData::controlTimes()
{
    scalarField controlTimes(totalNContSteps_+1, 0.0);

    const scalar& startTime = time_.startTime().value();

    forAll(controlTimes, tT)
    {
        controlTimes[tT] = startTime + tT*controlTime_.deltaT();
    }

    return controlTimes;
}

scalarField timeFluxData::averagingTimes()
{
    scalarField averagingTimes(totalNAvSteps_+1, 0.0);

    const scalar& startTime = time_.startTime().value();

    forAll(averagingTimes, tT)
    {
        averagingTimes[tT] = startTime + tT*averagingTime_.deltaT();
    }

    return averagingTimes;
}

scalarField timeFluxData::writeTimes()
{
    scalarField writeTimes(totalNWrSteps_+1, 0.0);

    const scalar& startTime = time_.startTime().value();

    forAll(writeTimes, tT)
    {
        writeTimes[tT] = startTime + tT*writeTime_.deltaT();
    }

    return writeTimes;
}




const timeInterval& timeFluxData::mdTimeInterval() const
{
    return mdTime_;
}

const timeInterval& timeFluxData::averageTimeInterval() const
{
    return averagingTime_;
}

const timeInterval& timeFluxData::controlTimeInterval() const
{
    return controlTime_;
}

const timeInterval& timeFluxData::writeTimeInterval() const
{
    return writeTime_;
}



timeInterval& timeFluxData::controlTimeInterval()
{
    return controlTime_;
}

timeInterval& timeFluxData::averageTimeInterval()
{
    return averagingTime_;
}

// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

//- Prefix increment
timeFluxData& timeFluxData::operator++()
{
    averagingTime_++;
    controlTime_++;
    writeTime_++;

    if(controlTime_.endTime())
    {
        controlTimeIndex_++;
    }

    if(averagingTime_.endTime())
    {
        averagingTimeIndex_++;
    }

    if(writeTime_.endTime())
    {
        writeTimeIndex_++;
    }

    return *this;
}


//- Postfix increment
timeFluxData& timeFluxData::operator++(int)
{
    return operator++();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
