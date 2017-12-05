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

#include "timeData.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Null Constructor
timeData::timeData
(
    Time& t
)
:
    time_(t),
    timeDict_(),
    writeInterval_(readScalar(t.controlDict().lookup("writeInterval"))),
    writeIntSteps_(label(writeInterval_/t.deltaT().value() + 0.5)),
    resetFieldsAtOutput_(true),
    resetIndex_(0),
    mdTime_(1),
    samplingTime_(),
    averagingTime_(),
//     calcPropTime_(),
    controlTime_(),
    nAvTimeSteps_("nAvTimeSteps_", dimless, 0.0),
    nControlSteps_(0.0),
    totalNSampSteps_(0),
    totalNAvSteps_(0),
    totalNContSteps_(0),
//     totalNCalcSteps_(0),
    controlTimeIndex_(0),
    averagingTimeIndex_(0),
//     calcTimeIndex_(0),
//     controlTimes_(),
//     calcTimes_()
//     averagingTimes_()

    decoupledFromWriteInterval_(false)
{}


//- Construct from Time and timeDict 
timeData::timeData
(
    Time& t,
    const dictionary& timeDict
)
:
    time_(t),
    timeDict_(timeDict),
    writeInterval_(readScalar(t.controlDict().lookup("writeInterval"))),
    writeIntSteps_(label(writeInterval_/t.deltaT().value() +0.5)),
    resetFieldsAtOutput_(true),
    resetIndex_(0),
    mdTime_(1),
    samplingTime_(readLabel(timeDict.lookup("nSamples"))),
    averagingTime_(readLabel(timeDict.lookup("nAverages"))),
//     calcPropTime_(readLabel(timeDict.lookup("nCalcProp"))),
    controlTime_(readLabel(timeDict.lookup("nControls"))),
    nAvTimeSteps_("nAvTimeSteps_", dimless, 0.0),
    nControlSteps_(0.0),
    totalNSampSteps_(0),
    totalNAvSteps_(0),
    totalNContSteps_(0),
//     totalNCalcSteps_(0),
    controlTimeIndex_(0),
//     calcTimeIndex_(0),
    averagingTimeIndex_(0),
    decoupledFromWriteInterval_(false)
//     controlTimes_(),
//     calcTimes_()
//     averagingTimes_()
{

    if (timeDict.found("decoupledFromWriteInterval"))
    {
        decoupledFromWriteInterval_ = Switch(timeDict.lookup("decoupledFromWriteInterval"));
    }


    setInitialData();

}

void timeData::setInitialData()
{
    Info << nl << "TimeData Statistics: " << endl;

    scalar deltaTMD = readScalar(time_.controlDict().lookup("deltaT"));

    mdTime_.deltaT() = deltaTMD;

    if (timeDict_.found("resetAtOutput"))
    {
        resetFieldsAtOutput_ = Switch(timeDict_.lookup("resetAtOutput"));
    }

    if(!decoupledFromWriteInterval_)
    {
        checkAndModifyTimeProperties();
    }
    else
    {
        averagingTime_.deltaT() = deltaTMD * averagingTime_.nSteps();
        controlTime_.deltaT() = deltaTMD * controlTime_.nSteps();
        samplingTime_.deltaT() = deltaTMD * samplingTime_.nSteps();
    }

    scalar endTime = time_.endTime().value();
    const scalar& startTime = time_.startTime().value();

    totalNAvSteps_ = label(((endTime - startTime) / averagingTime_.deltaT())+0.5);

    totalNContSteps_ = label(((endTime - startTime) / controlTime_.deltaT())+0.5);

//     totalNCalcSteps_ = label((endTime - startTime) / calcPropTime_.deltaT());

    totalNSampSteps_ = label(((endTime - startTime) / samplingTime_.deltaT())+0.5);

    Info << " total no. of sampling steps: " << totalNSampSteps_ << endl;
    Info << " total no. of averaging Steps: " << totalNAvSteps_ << endl;
    Info << " total no. of control Steps: " << totalNContSteps_ << endl;
//     Info << " total no. of calc Steps: " << totalNCalcSteps_ << endl;

    Info << nl << endl;

//     controlTimes_.setSize(totalNContSteps_+1, 0.0);
// //     calcTimes_.setSize(totalNCalcSteps_+1, 0.0);
//     averagingTimes_.setSize(totalNAvSteps_+1, 0.0);


//     forAll(controlTimes_, tT)
//     {
//         controlTimes_[tT] = startTime + tT*controlTime_.deltaT();
//     }

//     forAll(calcTimes_, tT)
//     {
//         calcTimes_[tT] = startTime + tT*calcPropTime_.deltaT();
//     }

//     forAll(averagingTimes_, tT)
//     {
//         averagingTimes_[tT] = startTime + tT*averagingTime_.deltaT();
//     }

    nAvTimeSteps_.value() = scalar(averagingTime_.nSteps());


    nControlSteps_ = averagingTime_.deltaT()/controlTime_.deltaT();


   //-- offsetting the controlling time index so that the time-interval finishes 
    //   one time-step ahead of the calcProp time.
    controlTime_.timeIndex()--;
}


void timeData::checkAndModifyTimeProperties()
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

    averagingTime_.deltaT() = deltaTMD * scalar(nAverages);

    Info << ", (final): " << nAverages
         << " time interval: " << averagingTime_.deltaT()
         << endl;


    // - 2. calcProp time
    //   for now we ensure that the calculation property interval is either 
    //   equal to the writing interval or else smaller than and a divisible factor
    //   of the writing interval

//     label& nCalcProp = calcPropTime_.nSteps();
// 
//     Info << "  nCalcProp steps (initial): " << nCalcProp;
// 
//     if(nCalcProp < 1)
//     {
//         nCalcProp = 1;
//         changedProperties = true;
//     }
//     else
//     {
//         if(nCalcProp > writeIntSteps_)
//         {
//             nCalcProp = writeIntSteps_;
//             changedProperties = true;
//         }
//         else
//         {
//             while((writeIntSteps_ % nCalcProp) != 0)
//             {
//                 nCalcProp ++;
//                 changedProperties = true;
//             }
//         }
//     }
// 
//     calcPropTime_.deltaT() = deltaTMD * scalar(nCalcProp);
// 
//     Info << ", (final): " << nCalcProp
//          << " time interval: " << calcPropTime_.deltaT()
//          << endl;
// 

   // - 3. control time

    label& nControls = controlTime_.nSteps();

    Info << " nControls (initial): " << nControls;

    if(nControls < 1)
    {
        nControls = 1;
        changedProperties = true;
    }
    else
    {
        if(nControls > writeIntSteps_)
        {
            nControls = writeIntSteps_;
            changedProperties = true;
        }

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


    // - 4. sampling time
    label& nSamples = samplingTime_.nSteps();

    Info << " nSamples (initial): " << nSamples;

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

    samplingTime_.deltaT() = deltaTMD * scalar(nSamples);

    Info << ", (final): " << nSamples
         << " time interval: " << samplingTime_.deltaT()
         << endl;


    if(changedProperties)
    {
        FatalErrorIn("timeData::timeData()")
            << "Time data members have been changed."
            << " Check and change them appropriately from the time dictionary" 
            << nl 
            << exit(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

timeData::~timeData()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void timeData::setTimeData(const dictionary& timeDict)
{
    const label nSamples(readLabel(timeDict.lookup("nSamples")));
    const label nAverages(readLabel(timeDict.lookup("nAverages")));
//     const label nCalcProp(readLabel(timeDict.lookup("nCalcProp")));
    const label nControls(readLabel(timeDict.lookup("nControls")));
 
    samplingTime_.nSteps() = nSamples;
    averagingTime_.nSteps() = nAverages;
//     calcPropTime_.nSteps() = nCalcProp;
    controlTime_.nSteps() = nControls;

    setInitialData();
}


const Time& timeData::time() const
{
    return time_;
}

Time& timeData::time()
{
    return time_;
}

const scalar& timeData::writeInterval() const
{
    return writeInterval_;
}

const label& timeData::writeIntervalSteps() const
{
    return writeIntSteps_;
}

const bool& timeData::samplingTime() const
{
    return samplingTime_.endTime();
}

const bool& timeData::averagingTime() const
{
    return averagingTime_.endTime();
}

const bool& timeData::controlTime() const
{
    return controlTime_.endTime();
}

// const bool& timeData::calcPropTime() const
// {
//     return calcPropTime_.endTime();
// }
// 


const label& timeData::nSamples() const
{
    return samplingTime_.nSteps();
}

const label& timeData::nControls() const
{
    return controlTime_.nSteps();
}

const label& timeData::nAverages() const
{
    return averagingTime_.nSteps();
}

const dimensionedScalar& timeData::nAvTimeSteps() const
{
    return nAvTimeSteps_;
}

scalar timeData::nAveragingTimeSteps()
{
    return scalar(nAvTimeSteps().value()*resetIndex_);
}

const scalar& timeData::nControlSteps() const
{
    return nControlSteps_;
}
// const label& timeData::nCalcProp() const
// {
//     return calcPropTime_.nSteps();
// }

const label& timeData::totalNSampSteps() const
{
    return totalNSampSteps_;
}


const label& timeData::totalNAvSteps() const
{
    return totalNAvSteps_;
}

const label& timeData::totalNContSteps() const
{
    return totalNContSteps_;
}

// const label& timeData::totalNCalcSteps() const
// {
//     return totalNCalcSteps_;
// }

const label& timeData::controlTimeIndex() const
{
    return controlTimeIndex_;
}

// const label& timeData::calcTimeIndex() const
// {
//     return calcTimeIndex_;
// }

const label& timeData::averagingTimeIndex() const
{
    return averagingTimeIndex_;
}


// const scalarField& timeData::controlTimes() const
// {
//     return controlTimes_;
// }



// const scalarField& timeData::calcTimes() const
// {
//     return calcTimes_;
// }

// const scalarField& timeData::averagingTimes() const
// {
//     return averagingTimes_;
// }


scalarField timeData::controlTimes()
{
    scalarField controlTimes(totalNContSteps_+1, 0.0);

    const scalar& startTime = time_.startTime().value();

    forAll(controlTimes, tT)
    {
        controlTimes[tT] = startTime + tT*controlTime_.deltaT();
    }

    return controlTimes;
}

scalarField timeData::averagingTimes()
{
    scalarField averagingTimes(totalNAvSteps_+1, 0.0);

    const scalar& startTime = time_.startTime().value();

    forAll(averagingTimes, tT)
    {
        averagingTimes[tT] = startTime + tT*averagingTime_.deltaT();
    }

    return averagingTimes;
}


const timeInterval& timeData::mdTimeInterval() const
{
    return mdTime_;
}

const timeInterval& timeData::sampleTimeInterval() const
{
    return samplingTime_;
}

const timeInterval& timeData::averageTimeInterval() const
{
    return averagingTime_;
}

/*const timeInterval& timeData::calcPropTimeInterval() const
{
    return calcPropTime_;
}*/
const timeInterval& timeData::controlTimeInterval() const
{
    return controlTime_;
}

timeInterval& timeData::controlTimeInterval()
{
    return controlTime_;
}

timeInterval& timeData::averageTimeInterval()
{
    return averagingTime_;
}


const bool& timeData::resetFieldsAtOutput() const
{
    return resetFieldsAtOutput_;
}

bool& timeData::resetFieldsAtOutput()
{
    return resetFieldsAtOutput_;
}

/*timeInterval& timeData::calcPropTimeInterval()
{
    return calcPropTime_;
}*/
// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

//- Prefix increment
timeData& timeData::operator++()
{
    samplingTime_++;
    averagingTime_++;
    controlTime_++;

//     mdTime_++;

    if(controlTime_.endTime())
    {
        controlTimeIndex_++;
    }

    if(averagingTime_.endTime())
    {

        if(resetFieldsAtOutput_)
        {
            resetIndex_ = 0;
        }

        resetIndex_++;

        averagingTimeIndex_++;
    }

//     if(calcPropTime_.endTime())
//     {
//         calcTimeIndex_++;
//     }

    return *this;
}


//- Postfix increment
timeData& timeData::operator++(int)
{
    return operator++();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
