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

#include "couplingTimeData.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Null Constructor
couplingTimeData::couplingTimeData
(
    Time& continuumTime,
    Time& molecularTime
)
:
    continuumTime_(continuumTime),
    molecularTime_(molecularTime),
    couplingTimeDict_(),
    writeInterval_(readScalar(continuumTime_.controlDict().lookup("writeInterval"))),
    couplingTimeInterval_(0.0),
    molecularTimeInt_(),
    continuumTimeInt_(),
    molecularTimeOn_(),
    molecularTimeOff_(),
    continuumSolution_(false),
    molecularSolution_(false),
    molecularToContinuumBCs_(false),
    continuumToMolecularBCs_(false),
    couplingOption_(),
    decoupledScheme_(false),
    molecularSwitch_(true),
    runSequential_(true)

{}


//- Construct from Time and timeDict 
couplingTimeData::couplingTimeData
(
    Time& continuumTime,
    Time& molecularTime,
    const dictionary& timeDict
)
:
    continuumTime_(continuumTime),
    molecularTime_(molecularTime),
    couplingTimeDict_(timeDict),
    writeInterval_(readScalar(continuumTime_.controlDict().lookup("writeInterval"))),
    couplingTimeInterval_(readScalar(couplingTimeDict_.lookup("couplingTimeInterval"))),
    molecularTimeInt_
    (
        label
        (
            couplingTimeInterval_/
            readScalar(molecularTime.controlDict().lookup("deltaT"))
        ),
        couplingTimeInterval_
    ),
    continuumTimeInt_
    (
        label
        (
            couplingTimeInterval_/
            readScalar(continuumTime.controlDict().lookup("deltaT"))
        ),
        couplingTimeInterval_
    ),
    molecularTimeOn_
    (
        label
        (
            couplingTimeInterval_/
            readScalar(molecularTime.controlDict().lookup("deltaT"))
        ),
        couplingTimeInterval_
    ),
    molecularTimeOff_
    (
        label
        (
            couplingTimeInterval_/
            readScalar(molecularTime.controlDict().lookup("deltaT"))
        ),
        couplingTimeInterval_    
    ),
    continuumSolution_(false),
    molecularSolution_(false),
    couplingOption_(couplingTimeDict_.lookup("couplingOption")),
    decoupledScheme_(false),
    molecularSwitch_(true),
    runSequential_(true)
{

    scalar deltaTC = readScalar(continuumTime_.controlDict().lookup("deltaT"));
    scalar deltaTMD = readScalar(molecularTime_.controlDict().lookup("deltaT"));

    //- checks

    //- coupling option: sequential, synchronous or decoupled

    if(couplingOption_ == "sequential")
    {
        continuumSolution_ = true;
        molecularSolution_ = false;
    }
    else if(couplingOption_ == "synchronous")
    {
        runSequential_ = Switch(couplingTimeDict_.lookup("runSequential"));

        if(runSequential_) // continuum "block" first
        {
            continuumSolution_ = true;
            molecularSolution_ = false;
        }
        else // continuum and molecular MD time-steps proceed
        {
            continuumSolution_ = true;
            molecularSolution_ = true;
        }
    }
    else if(couplingOption_ == "decoupled")
    {
        decoupledScheme_ = true;
        molecularSwitch_ = false;
        continuumSolution_ = true;
        molecularSolution_ = false;

        scalar molecularOnInterval = readScalar(couplingTimeDict_.lookup("molecularOnInterval"));

        if(molecularOnInterval > couplingTimeInterval_)
        {
            FatalErrorIn("couplingTimeData::couplingTimeData()")
                << "molecularOnInterval: " <<  molecularOnInterval 
                << "should be less than couplingTimeInterval: " << couplingTimeInterval_
                << nl << continuumTime_.systemPath()/"couplingDict"
                << exit(FatalError);
        }

        if(molecularOnInterval <= 0.0)
        {
            FatalErrorIn("couplingTimeData::couplingTimeData()")
                << "molecularOnInterval: " <<  molecularOnInterval 
                << "should be greater than 0. "
                << nl << continuumTime_.systemPath()/"couplingDict"
                << exit(FatalError);
        }

        scalar molecularOffInterval = couplingTimeInterval_ - molecularOnInterval;

        if(molecularOffInterval == 0.0)
        {
            molecularSwitch_ = true;
            decoupledScheme_ = false;
            couplingOption_ = "sequential";

            Info << " changing coupling option from: decoupled to: sequential" << endl;
        }

        label nStepsOn = label(molecularOnInterval/deltaTMD);

        molecularTimeOn_.deltaT() = scalar(nStepsOn)*deltaTMD;
        molecularTimeOn_.nSteps() = nStepsOn;

        label nStepsOff = label(molecularOffInterval/deltaTMD);

        molecularTimeOff_.deltaT() = scalar(nStepsOff)*deltaTMD;
        molecularTimeOff_.nSteps() = nStepsOff;
    }
    else 
    {
        FatalErrorIn("couplingTimeData::couplingTimeData()")
            << "Coupling options available: sequential or sequential in " 
            << nl << continuumTime_.systemPath()/"couplingDict"
            << exit(FatalError);
    }

    //- write interval
    scalar molecularWriteInterval = readScalar(molecularTime_.controlDict().lookup("writeInterval"));

    if(molecularWriteInterval != writeInterval_)
    {
        FatalErrorIn("couplingTimeData::couplingTimeData()")
            << "Write intervals must be the same for both solvers: " << nl
            << "Molecular write interval is: " << molecularWriteInterval << nl
            << "Continuum write interval is: " << writeInterval_ << nl
            << exit(FatalError);
    }

    Info << nl << "Time-Coupling Statistics: " << endl;

    Info << " continuum time-step = " << deltaTC << endl;
    Info << " molecular time-step = " << deltaTMD << endl;
    Info << " time-step disparity ratio = " << couplingRatio() << endl;
    Info << " coupling time-interval = " << couplingTimeInterval_ << endl;
    Info << " coupling option: " << couplingOption_ << nl << endl;

    if(couplingOption_ == "decoupled")
    {
        Info << " molecular-time interval on = " << molecularTimeOn_.deltaT() << endl;
        Info << " molecular-time interval off = " << molecularTimeOff_.deltaT() << endl;

        Info << " decoupled ratio = " 
             << couplingTimeInterval_/molecularTimeOn_.deltaT() 
             << nl << endl;
    }

    scalar continuumError = 
    mag
    (
        scalar(continuumTimeInt_.nSteps()) - (couplingTimeInterval_/deltaTC)
    );

    scalar molecularError = 
    mag
    (
        scalar(molecularTimeInt_.nSteps()) - (couplingTimeInterval_/deltaTMD)
    );

    if(continuumError > SMALL)
    {
        FatalErrorIn("couplingTimeData::couplingTimeData()")
            << "The ratio: (coupling time interval / continuum time), must be an integer: "
            << couplingTimeInterval_ << " / " << deltaTC << " = "
            << couplingTimeInterval_/deltaTC << nl
            << exit(FatalError);
    }

    if(molecularError > SMALL)
    {
        FatalErrorIn("couplingTimeData::couplingTimeData()")
            << "The ratio: (coupling time interval / molecular time), must be an integer: "
            << couplingTimeInterval_ << " / " << deltaTMD << " = "
            << couplingTimeInterval_/deltaTMD << nl
            << exit(FatalError);
    }

    if(couplingTimeInterval_ < deltaTC)
    {
        FatalErrorIn("couplingTimeData::couplingTimeData()")
            << "Choose the coupling time interval greater or equal to the continuum time-interval." << nl
            << exit(FatalError);
    }

    if(couplingTimeInterval_ < deltaTMD)
    {
        FatalErrorIn("couplingTimeData::couplingTimeData()")
            << "Choose the coupling time interval greater or equal to the molecular time-interval." << nl
            << exit(FatalError);
    }

    if(couplingTimeInterval_ > writeInterval_)
    {
        FatalErrorIn("couplingTimeData::couplingTimeData()")
            << "Choose the coupling time interval less or equal to the write interval"
            << nl
            << exit(FatalError);
    }


}

// void couplingTimeData::setInitialData()
// {
//     Info << nl << "TimeData Statistics: " << endl;
// 
//     scalar deltaTMD = readScalar(time_.controlDict().lookup("deltaT"));
// 
//     mdTime_.deltaT() = deltaTMD;
// 
//     if(timeMeasOption_ == "write")
//     {
//         samplingTime_.nSteps() = 1;
//         averagingTime_.nSteps() = writeIntSteps_;
//     }
//     else if(timeMeasOption_ == "general")
//     {
//         const label nSamples = readLabel(timeDict_.lookup("nSamples"));
//         const label nAverages = readLabel(timeDict_.lookup("nAverages"));
// 
//         samplingTime_.nSteps() = nSamples;
//         averagingTime_.nSteps() = nAverages;
// 
//         checkAndModifyTimeProperties();
//     }
//     else
//     {
//         FatalErrorIn("couplingTimeData::setInitialData()")
//             << "timeOption: \"" << timeMeasOption_
//             << "\" is not one of the available options." << nl 
//             << exit(FatalError);
//     }
// 
//     samplingTime_.deltaT() = deltaTMD * scalar(samplingTime_.nSteps());
//     averagingTime_.deltaT() = deltaTMD * scalar(averagingTime_.nSteps());
// 
//     Info << " nSamples: " << samplingTime_.nSteps()
//          << ", time interval: " << samplingTime_.deltaT()
//          << endl;
// 
//     Info << " nAverages: " << averagingTime_.nSteps()
//          << ", time interval: " << averagingTime_.deltaT()
//          << endl;
// 
// 
//     const scalar& endTime = time_.endTime().value();
//     const scalar& startTime = time_.startTime().value();
// 
//     totalNAvSteps_ = label((endTime - startTime) / averagingTime_.deltaT());
// 
// //     totalNCalcSteps_ = label((endTime - startTime) / calcPropTime_.deltaT());
// 
//     totalNSampSteps_ = label((endTime - startTime) / samplingTime_.deltaT());
// 
//     Info << " total no. of sampling steps: " << totalNSampSteps_ << endl;
//     Info << " total no. of averaging Steps: " << totalNAvSteps_ << endl;
// //     Info << " total no. of calc Steps: " << totalNCalcSteps_ << endl;
// 
//     Info << nl << endl;
// 
//     averagingTimes_.setSize(totalNAvSteps_+1, 0.0);
// 
//     forAll(averagingTimes_, tT)
//     {
//         averagingTimes_[tT] = startTime + tT*averagingTime_.deltaT();
//     }
// 
//     nAvTimeSteps_.value() = scalar(averagingTime_.nSteps());
// 
//     samplingTimes_.setSize(totalNSampSteps_+1, 0.0);
// 
//     forAll(samplingTimes_, tT)
//     {
//         samplingTimes_[tT] = startTime + tT*samplingTime_.deltaT();
//     }
// 
// }


// void couplingTimeData::checkAndModifyTimeProperties()
// {
//     //- checking 
// 
//     bool changedProperties = false;
// 
//     const label nAveragesOriginal = averagingTime_.nSteps();
//     const label nSamplesOriginal = samplingTime_.nSteps();
// 
//     // - 1. averaging time
// 
//     label& nAverages = averagingTime_.nSteps();
// 
//     if(nAverages < 1)
//     {
//         nAverages = 1;
//         changedProperties = true;
//     }
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
// 
//     // - 2. sampling time
//     label& nSamples = samplingTime_.nSteps();
// 
//     if(nSamples < 1)
//     {
//         nSamples = 1;
//         changedProperties = true;
//     }
//     else
//     {
//         if(nSamples > nAverages)
//         {
//             nSamples = nAverages;
//             changedProperties = true;
//         }
//         else
//         {
//             while((nAverages % nSamples) != 0)
//             {
//                 nSamples--;
//                 changedProperties = true;
//             }
//         }
//     }
// 
//     if(changedProperties)
//     {
//         Info << "initial nSamples: " << nSamplesOriginal
//              << ", modified nSamples: " << samplingTime_.nSteps() << endl;
// 
//         Info << "initial nAverages: " << nAveragesOriginal
//              << ", modified nAverages: " << averagingTime_.nSteps() << endl;
// 
//         FatalErrorIn("couplingTimeData::couplingTimeData()")
//             << "Time properties are inconsistent."
//             << " Check and change them appropriately from the time dictionary" 
//             << nl 
//             << exit(FatalError);
//     }
// }


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

couplingTimeData::~couplingTimeData()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// void couplingTimeData::setTimeData(const dictionary& timeDict)
// {
//     const label nSamples(readLabel(timeDict.lookup("nSamples")));
//     const label nAverages(readLabel(timeDict.lookup("nAverages")));
// //     const label nCalcProp(readLabel(timeDict.lookup("nCalcProp")));
//     
//  
//     samplingTime_.nSteps() = nSamples;
//     averagingTime_.nSteps() = nAverages;
// //     calcPropTime_.nSteps() = nCalcProp;
//     
// 
//     setInitialData();
// }


const Time& couplingTimeData::continuumTime() const
{
    return continuumTime_;
}

Time& couplingTimeData::continuumTime()
{
    return continuumTime_;
}

const Time& couplingTimeData::molecularTime() const
{
    return molecularTime_;
}

Time& couplingTimeData::molecularTime()
{
    return molecularTime_;
}


const scalar& couplingTimeData::writeInterval() const
{
    return writeInterval_;
}


const bool& couplingTimeData::continuumSolution() const
{
    return continuumSolution_;
}

bool& couplingTimeData::continuumSolution()
{
    return continuumSolution_;
}

const bool& couplingTimeData::molecularSolution() const
{
    return molecularSolution_;
}

bool& couplingTimeData::molecularSolution()
{
    return molecularSolution_;
}

// const label& couplingTimeData::writeIntervalSteps() const
// {
//     return writeIntSteps_;
// }

// const bool& couplingTimeData::samplingTime() const
// {
//     return samplingTime_.endTime();
// }
// 
// const bool& couplingTimeData::averagingTime() const
// {
//     return averagingTime_.endTime();
// }

// const bool& couplingTimeData::calcPropTime() const
// {
//     return calcPropTime_.endTime();
// }
// 


// const label& couplingTimeData::nSamples() const
// {
//     return samplingTime_.nSteps();
// }
// 
// 
// const label& couplingTimeData::nAverages() const
// {
//     return averagingTime_.nSteps();
// }
// 
// const dimensionedScalar& couplingTimeData::nAvTimeSteps() const
// {
//     return nAvTimeSteps_;
// }
// 
// // const label& couplingTimeData::nCalcProp() const
// // {
// //     return calcPropTime_.nSteps();
// // }
// 
// const label& couplingTimeData::totalNSampSteps() const
// {
//     return totalNSampSteps_;
// }
// 
// const label& couplingTimeData::totalNAvSteps() const
// {
//     return totalNAvSteps_;
// }
// 
// const label& couplingTimeData::averagingTimeIndex() const
// {
//     return averagingTimeIndex_;
// }
// 
// const label& couplingTimeData::samplingTimeIndex() const
// {
//     return samplingTimeIndex_;
// }
// 
// // const label& couplingTimeData::totalNCalcSteps() const
// // {
// //     return totalNCalcSteps_;
// // }
// 
// // const label& couplingTimeData::calcTimeIndex() const
// // {
// //     return calcTimeIndex_;
// // }
// 
// const scalarField& couplingTimeData::averagingTimes() const
// {
//     return averagingTimes_;
// }
// const scalarField& couplingTimeData::samplingTimes() const
// {
//     return samplingTimes_;
// }
// 
// const timeInterval& couplingTimeData::mdTimeInterval() const
// {
//     return mdTime_;
// }
// 
// const timeInterval& couplingTimeData::sampleTimeInterval() const
// {
//     return samplingTime_;
// }
// 
// const timeInterval& couplingTimeData::averageTimeInterval() const
// {
//     return averagingTime_;
// }


// const timeInterval& couplingTimeData::calcPropTimeInterval() const
// {
//     return calcPropTime_;
// }

// timeInterval& couplingTimeData::calcPropTimeInterval()
// {
//     return calcPropTime_;
// }

scalar couplingTimeData::couplingRatio()
{
    scalar deltaTC = readScalar(continuumTime_.controlDict().lookup("deltaT"));
    scalar deltaTMD = readScalar(molecularTime_.controlDict().lookup("deltaT"));

    return deltaTC/deltaTMD;
}

const timeInterval& couplingTimeData::continuumTimeInterval() const
{
    return continuumTimeInt_;
}

const word& couplingTimeData::couplingOption() const
{
    return couplingOption_;
}

const timeInterval& couplingTimeData::molecularTimeInterval() const
{
    return molecularTimeInt_;
}

const bool& couplingTimeData::continuumToMolecularBCs() const
{
    return continuumToMolecularBCs_;
}

const bool& couplingTimeData::molecularToContinuumBCs() const
{
    return molecularToContinuumBCs_;    
}

const bool& couplingTimeData::molecularSwitch() const
{
    return molecularSwitch_;    
}


//- place this function right at the end of the hybrid solver while loop
void couplingTimeData::solversToRun()
{
    if(couplingOption_ == "sequential")
    {
        //-switch solvers 

        if(molecularTimeInt_.endTime())
        {
            molecularSolution_ = false;
            continuumSolution_ = true;

            //- temporary
            molecularTimeInt_.endTime() = false;
        }
    
        if(continuumTimeInt_.endTime())
        {
            continuumSolution_ = false;
            molecularSolution_ = true;

            continuumTimeInt_.endTime() = false;
        }
    }

    if(couplingOption_ == "synchronous")
    {
        if
        (
            (!molecularTimeInt_.endTime()) && 
            (continuumTimeInt_.endTime())
        )
        {
            continuumSolution_ = false;
            molecularSolution_ = true;
        }

        // - just in case 
        if
        (
            (molecularTimeInt_.endTime()) && 
            (!continuumTimeInt_.endTime())
        )
        {
            continuumSolution_ = true;
            molecularSolution_ = false;
        }

        if
        (
            (molecularTimeInt_.endTime()) && 
            (continuumTimeInt_.endTime())
        )
        {
            if(runSequential_) // continuum "block" first
            {
                continuumSolution_ = true;
                molecularSolution_ = false;
            }
            else // continuum and molecular MD time-steps proceed
            {
                continuumSolution_ = true;
                molecularSolution_ = true;
            }

            continuumTimeInt_.endTime() = false;
            molecularTimeInt_.endTime() = false;
        }
    }


    if(couplingOption_ == "decoupled")
    {
        if(continuumTimeInt_.endTime())
        {
            continuumSolution_ = false;
            molecularSolution_ = true;
            molecularSwitch_ = false;
            continuumTimeInt_.endTime() = false;
        }

        if(molecularTimeInt_.endTime())
        {
            molecularSolution_ = false;
            continuumSolution_ = true;
            molecularSwitch_ = false;

            //- temporary
            molecularTimeInt_.endTime() = false;
        }


        if(molecularTimeOff_.endTime())
        {
            molecularSwitch_ = true;

            molecularTimeOff_.endTime() = false;
        }
    }
}


void couplingTimeData::boundaryConditions()
{
    if(couplingOption_ == "sequential")
    {
        //- boundary conditions tested here
        if(molecularTimeInt_.endTime())
        {
            molecularToContinuumBCs_ = true;
        }
        else
        {
            molecularToContinuumBCs_ = false;
        }
    
        if(continuumTimeInt_.endTime())
        {
            continuumToMolecularBCs_ = true;
        }
        else
        {
            continuumToMolecularBCs_ = false;
        }
    }

    if(couplingOption_ == "synchronous")
    {
        //- boundary conditions tested here
        if
        (
            (molecularTimeInt_.endTime()) &&
            (continuumTimeInt_.endTime()) 
        )
        {
            molecularToContinuumBCs_ = true;
            continuumToMolecularBCs_ = true;
        }
        else
        {
            molecularToContinuumBCs_ = false;
            continuumToMolecularBCs_ = false;
        }
    }

    if(couplingOption_ == "decoupled") // same as a sequential scheme
    {
        //- boundary conditions tested here
        if(molecularTimeInt_.endTime())
        {
            molecularToContinuumBCs_ = true;
        }
        else
        {
            molecularToContinuumBCs_ = false;
        }
    
        if(continuumTimeInt_.endTime())
        {
            continuumToMolecularBCs_ = true;
        }
        else
        {
            continuumToMolecularBCs_ = false;
        }
    }
}

// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

//- Prefix increment
couplingTimeData& couplingTimeData::operator++()
{
    if(molecularSolution_)
    {
        molecularTimeInt_++;

        if(decoupledScheme_ && !molecularSwitch_)
        {
            molecularTimeOff_++;
        }

        if(decoupledScheme_ && molecularSwitch_)
        {
            molecularTimeOn_++;
        }
    }

    if(continuumSolution_)
    {
        continuumTimeInt_++;
    }


    boundaryConditions();

    return *this;
}


//- Postfix increment
couplingTimeData& couplingTimeData::operator++(int)
{
    return operator++();
}



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
