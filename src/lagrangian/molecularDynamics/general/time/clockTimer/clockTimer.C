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

#include "clockTimer.H"
#include "writeTimeData.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

//-  Constructor
clockTimer::clockTimer(Time& t, word fieldName, bool write)
:
    time_(t),
    fieldName_(fieldName),
    instTimeIndex_(0.0),
    timeIndex_(0.0),
    totalDuration_(0.0),
    duration_(0.0),
    instantDuration_(0.0),
    writeOut_(write),
    timePath_()
{
    if(writeOut_)
    {
        // directory: case/boundaries
        timePath_ = time_.path()/"timings";

        if(isDir(timePath_))
        {
            rmDir(timePath_);
        }

        if(Pstream::master())
        {
            mkDir(timePath_);
        }
    }



}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

clockTimer::~clockTimer()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
// void clockTimer::startClock()
// {
//     scalar newTime = time_.elapsedCpuTime();
// 
//     if(timeIndex_ == 0.0)
//     {
//         startTimeBefore_ = newTime;
// 
//         Info << "startTimeBefore: " << startTimeBefore_ << endl;
//     }
//     else
//     {
//         innerIncrement_ += newTime - startTimeAfter_;
//     }
// 
// //     Info << "time before: " << lastTime_ << endl;
// }
// 
// void clockTimer::stopClock()
// {
//     scalar newTime = time_.elapsedCpuTime();
// 
//     if(timeIndex_ == 0.0)
//     {
//         startTimeAfter_ = newTime;
// 
//         Info << "startTimeAfter: " << startTimeAfter_ << endl;
//     }
//     else
//     {
//         outerIncrement_ = newTime - startTimeBefore_;
//     }
// 
//     timeIndex_ += 1.0;
//     
//     Info << "outerIncrement: " << outerIncrement_ << endl;
//     Info << "innerIncrement: " << innerIncrement_ << endl;
// 
//     Info << "average timing: " << averageTime() << endl;
// }

void clockTimer::startClock()
{
    gettimeofday(&lastTime_, NULL);
//     clock_gettime(CLOCK_MONOTONIC, &lastTime_);
}

void clockTimer::stopClock()
{
    timeval endTime;
//     timespec endTime;
    scalar seconds = 0.0;
    scalar useconds = 0.0;

//     clock_gettime(CLOCK_MONOTONIC, &endTime);
    gettimeofday(&endTime, NULL);

    seconds  = endTime.tv_sec  - lastTime_.tv_sec;
//     nseconds = endTime.tv_nsec - lastTime_.tv_nsec;
    useconds = endTime.tv_usec - lastTime_.tv_usec;

    scalar duration = seconds + (useconds*1e-6);

    duration_ = duration;
    instantDuration_ += duration;
    totalDuration_ += duration;

    instTimeIndex_ += 1.0;
    timeIndex_ += 1.0;

    
    Info<< "Duration: " << fieldName_ << ", inst. = " << duration
        << " s   av. write int. = " << averageTimeWriteInterval()
        << " s   av. sim. = " << averageTime()
        << " s   tot. = " << totalDuration_ << " s"
        << endl;

//     elapsedMicroCpuField_[index_] = duration_;
//     elapsedCpuField_[index_] = time_.elapsedCpuTime();



    write();
}

// scalar clockTimer::averageTime()
// {
//     return (outerIncrement_-innerIncrement_)/timeIndex_;
// }
    
    
scalar clockTimer::averageTimeWriteInterval()
{
    if(instTimeIndex_ > 0)
    {
        return instantDuration_/instTimeIndex_;
    }
    else
    {
        return 0.0;
    }
}

scalar clockTimer::averageTime()
{
    if(timeIndex_ > 0)
    {
        return totalDuration_/timeIndex_;
    }
    else
    {
        return 0.0;
    }
}


const scalar& clockTimer::instantDuration() const
{
    return duration_;
}

    
const scalar& clockTimer::totalDuration() const
{
    return totalDuration_;
}

void clockTimer::write()
{
    if(time_.outputTime())
    {

        if(Pstream::master() && writeOut_)
        {
            scalarField timeField(1, time_.timeOutputValue());
            scalarField instantDuration(1, averageTimeWriteInterval());
            scalarField cumulDuration(1, averageTime());

//             fileName casePath(time_.path());

//             writeTimeData
//             (
//                 timePath_,
//                 "cpuTimeProcess_"+fieldName_+"_instant.xy",
//                 timeField,
//                 instantTimeDuration_,
//                 true
//             );

            writeTimeData
            (
                 timePath_,
                 "cpuTimeProcess_"+fieldName_+"_instant.xy",
                 timeField,
                 instantDuration,
                 true
             );
            
            writeTimeData
            (
                timePath_,
                "cpuTimeProcess_"+fieldName_+"_average.xy",
                timeField,
                cumulDuration,
                true
            );

        }
        
        // reset
        instTimeIndex_ = 0.0;
        instantDuration_ = 0.0;
    }
}

// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
