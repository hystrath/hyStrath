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

#include "decoupledTimeData.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Null Constructor
decoupledTimeData::decoupledTimeData
(
    Time& continuumTime,
    Time& molecularTime
)
:
    continuumTime_(continuumTime),
    molecularTime_(molecularTime),
    couplingTimeDict_(),

    contWriteInterval_(readScalar(continuumTime_.controlDict().lookup("writeInterval"))),
    molWriteInterval_(readScalar(molecularTime_.controlDict().lookup("writeInterval"))),

    contCouplingTimeInterval_(0.0),
    molCouplingTimeInterval_(0.0),

    molecularTimeInt_(),
    continuumTimeInt_(),

    continuumSolution_(false),
    molecularSolution_(false),

    molecularToContinuumBCs_(false),
    continuumToMolecularBCs_(false),

    couplingOption_()/*
    decoupledScheme_(false),
    molecularSwitch_(true),
    runSequential_(true)*/

{}


//- Construct from Time and timeDict 
decoupledTimeData::decoupledTimeData
(
    Time& continuumTime,
    Time& molecularTime,
    const dictionary& timeDict
)
:
    continuumTime_(continuumTime),
    molecularTime_(molecularTime),
    couplingTimeDict_(timeDict),
    contWriteInterval_(readScalar(continuumTime_.controlDict().lookup("writeInterval"))),
    molWriteInterval_(readScalar(molecularTime_.controlDict().lookup("writeInterval"))),
    contCouplingTimeInterval_(readScalar(couplingTimeDict_.lookup("continuumCouplingTimeInterval"))),
    molCouplingTimeInterval_(readScalar(couplingTimeDict_.lookup("molecularCouplingTimeInterval"))),
    molecularTimeInt_
    (
        label
        (
            molCouplingTimeInterval_/
            readScalar(molecularTime.controlDict().lookup("deltaT"))
        ),
        molCouplingTimeInterval_
    ),
    continuumTimeInt_
    (
        label
        (
            contCouplingTimeInterval_/
            readScalar(continuumTime.controlDict().lookup("deltaT"))
        ),
        contCouplingTimeInterval_
    ),
    continuumSolution_(false),
    molecularSolution_(false),
    molecularToContinuumBCs_(false),
    continuumToMolecularBCs_(true),
    couplingOption_(couplingTimeDict_.lookup("couplingOption"))
//     decoupledScheme_(false),
//     molecularSwitch_(true),
//     runSequential_(true)
{

    scalar deltaTC = readScalar(continuumTime_.controlDict().lookup("deltaT"));
    scalar deltaTMD = readScalar(molecularTime_.controlDict().lookup("deltaT"));

    //- checks
    //- coupling option: sequential, synchronous or decoupled

    if(couplingOption_ == "sequential")
    {
        continuumSolution_ = false;
        molecularSolution_ = true;
    }
    else 
    {
        FatalErrorIn("decoupledTimeData::decoupledTimeData()")
            << "Coupling options available: sequential" 
            << nl << " in : " 
            << continuumTime_.systemPath()/"couplingDict"
            << exit(FatalError);
    }


//     else if(couplingOption_ == "decoupled")
//     {
//         decoupledScheme_ = true;
//         molecularSwitch_ = false;
//         continuumSolution_ = true;
//         molecularSolution_ = false;
// 
//         scalar molecularOnInterval = readScalar(couplingTimeDict_.lookup("molecularOnInterval"));
// 
//         if(molecularOnInterval > couplingTimeInterval_)
//         {
//             FatalErrorIn("decoupledTimeData::decoupledTimeData()")
//                 << "molecularOnInterval: " <<  molecularOnInterval 
//                 << "should be less than couplingTimeInterval: " << couplingTimeInterval_
//                 << nl << continuumTime_.systemPath()/"couplingDict"
//                 << exit(FatalError);
//         }
// 
//         if(molecularOnInterval <= 0.0)
//         {
//             FatalErrorIn("decoupledTimeData::decoupledTimeData()")
//                 << "molecularOnInterval: " <<  molecularOnInterval 
//                 << "should be greater than 0. "
//                 << nl << continuumTime_.systemPath()/"couplingDict"
//                 << exit(FatalError);
//         }
// 
//         scalar molecularOffInterval = couplingTimeInterval_ - molecularOnInterval;
// 
//         if(molecularOffInterval == 0.0)
//         {
//             molecularSwitch_ = true;
//             decoupledScheme_ = false;
//             couplingOption_ = "sequential";
// 
//             Info << " changing coupling option from: decoupled to: sequential" << endl;
//         }
// 
//         label nStepsOn = label(molecularOnInterval/deltaTMD);
// 
//         molecularTimeOn_.deltaT() = scalar(nStepsOn)*deltaTMD;
//         molecularTimeOn_.nSteps() = nStepsOn;
// 
//         label nStepsOff = label(molecularOffInterval/deltaTMD);
// 
//         molecularTimeOff_.deltaT() = scalar(nStepsOff)*deltaTMD;
//         molecularTimeOff_.nSteps() = nStepsOff;
//     }
//     else 
//     {
//         FatalErrorIn("decoupledTimeData::decoupledTimeData()")
//             << "Coupling options available: sequential or sequential in " 
//             << nl << continuumTime_.systemPath()/"couplingDict"
//             << exit(FatalError);
//     }

    //- write interval
//     scalar molecularWriteInterval = readScalar(molecularTime_.controlDict().lookup("writeInterval"));
// 
//     if(molecularWriteInterval != writeInterval_)
//     {
//         FatalErrorIn("decoupledTimeData::decoupledTimeData()")
//             << "Write intervals must be the same for both solvers: " << nl
//             << "Molecular write interval is: " << molecularWriteInterval << nl
//             << "Continuum write interval is: " << writeInterval_ << nl
//             << exit(FatalError);
//     }

    Info << nl << "Hybrid molecular-continuum time-coupling statistics: " << nl << endl;

    Info << "  continuum time-step = " << deltaTC << endl;
    Info << "  molecular time-step = " << deltaTMD << endl;
    Info << "  time-step disparity ratio = " << couplingRatio() << endl;
    Info << "  continuum coupling time-interval = " << contCouplingTimeInterval_ << endl;
    Info << "  molecular coupling time-interval = " << molCouplingTimeInterval_ << endl;
    Info << "  decoupling ratio = " << contCouplingTimeInterval_/molCouplingTimeInterval_ << endl;
    Info << "  coupling option: " << couplingOption_ << nl << endl;



    //- Checking for consistency -- ensuring that both molecular and continuum time-lines 
    // end at the same time.


    const scalar& endTimeCont = continuumTime_.endTime().value();
    const scalar& startTimeCont = continuumTime_.startTime().value();

    const scalar& endTimeMol = molecularTime_.endTime().value();
    const scalar& startTimeMol = molecularTime_.startTime().value();

    label nCouplingStepsCont = label((endTimeCont - startTimeCont)/contCouplingTimeInterval_);
    label nCouplingStepsMol = label((endTimeMol - startTimeMol)/molCouplingTimeInterval_);


    if(nCouplingStepsCont != nCouplingStepsMol)
    {
        scalar endTimeMolMod = startTimeMol + (molCouplingTimeInterval_/contCouplingTimeInterval_)*(endTimeCont - startTimeCont);

        FatalErrorIn("decoupledTimeData::decoupledTimeData()")
            << "Inconsistency in time-scheme. "
            << "The total number of coupling steps between the continuum: " 
            << nCouplingStepsCont << ", and the molecular: " 
            << nCouplingStepsMol << ", should be equal."
            << " Check start/end times, and coupling time-intervals. "
            << " Example: change molecular end-time from " 
            << endTimeMol << " to " << endTimeMolMod
            << exit(FatalError);
    }



//     if(couplingOption_ == "decoupled")
//     {
//         Info << " molecular-time interval on = " << molecularTimeOn_.deltaT() << endl;
//         Info << " molecular-time interval off = " << molecularTimeOff_.deltaT() << endl;
// 
//         Info << " decoupled ratio = " 
//              << couplingTimeInterval_/molecularTimeOn_.deltaT() 
//              << nl << endl;
//     }

//     scalar continuumError = 
//     mag
//     (
//         scalar(continuumTimeInt_.nSteps()) - (couplingTimeInterval_/deltaTC)
//     );
// 
//     scalar molecularError = 
//     mag
//     (
//         scalar(molecularTimeInt_.nSteps()) - (couplingTimeInterval_/deltaTMD)
//     );
// 
//     if(continuumError > SMALL)
//     {
//         FatalErrorIn("decoupledTimeData::decoupledTimeData()")
//             << "The ratio: (coupling time interval / continuum time), must be an integer: "
//             << couplingTimeInterval_ << " / " << deltaTC << " = "
//             << couplingTimeInterval_/deltaTC << nl
//             << exit(FatalError);
//     }
// 
//     if(molecularError > SMALL)
//     {
//         FatalErrorIn("decoupledTimeData::decoupledTimeData()")
//             << "The ratio: (coupling time interval / molecular time), must be an integer: "
//             << couplingTimeInterval_ << " / " << deltaTMD << " = "
//             << couplingTimeInterval_/deltaTMD << nl
//             << exit(FatalError);
//     }

    if(contCouplingTimeInterval_ < deltaTC)
    {
        FatalErrorIn("decoupledTimeData::decoupledTimeData()")
            << "Choose the continuum coupling time interval: " 
            << contCouplingTimeInterval_
            << " greater or equal to the continuum time-step: " 
            << deltaTC
            << nl
            << exit(FatalError);
    }

    if(molCouplingTimeInterval_ < deltaTMD)
    {
        FatalErrorIn("decoupledTimeData::decoupledTimeData()")
            << "Choose the molecular coupling time interval" 
            << " greater or equal to the molecular time-interval."
            << nl
            << exit(FatalError);
    }

//     if(contCouplingTimeInterval_ > contWriteInterval_)
//     {
//         FatalErrorIn("decoupledTimeData::decoupledTimeData()")
//             << "Choose the continuum coupling time interval"
//             << " less or equal to the continuum, write interval"
//             << nl
//             << exit(FatalError);
//     }
// 
//     if(molCouplingTimeInterval_ > molWriteInterval_)
//     {
//         FatalErrorIn("decoupledTimeData::decoupledTimeData()")
//             << "Choose the molecular coupling time interval" 
//             << " less or equal to the molecular write interval"
//             << nl
//             << exit(FatalError);
//     }

}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

decoupledTimeData::~decoupledTimeData()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


const Time& decoupledTimeData::continuumTime() const
{
    return continuumTime_;
}

Time& decoupledTimeData::continuumTime()
{
    return continuumTime_;
}

const Time& decoupledTimeData::molecularTime() const
{
    return molecularTime_;
}

Time& decoupledTimeData::molecularTime()
{
    return molecularTime_;
}

const scalar& decoupledTimeData::contWriteInterval() const
{
    return contWriteInterval_;
}

const scalar& decoupledTimeData::molWriteInterval() const
{
    return molWriteInterval_;
}

const bool& decoupledTimeData::continuumSolution() const
{
    return continuumSolution_;
}

bool& decoupledTimeData::continuumSolution()
{
    return continuumSolution_;
}

const bool& decoupledTimeData::molecularSolution() const
{
    return molecularSolution_;
}

bool& decoupledTimeData::molecularSolution()
{
    return molecularSolution_;
}

scalar decoupledTimeData::couplingRatio()
{
    scalar deltaTC = readScalar(continuumTime_.controlDict().lookup("deltaT"));
    scalar deltaTMD = readScalar(molecularTime_.controlDict().lookup("deltaT"));

    return deltaTC/deltaTMD;
}

const timeInterval& decoupledTimeData::continuumTimeInterval() const
{
    return continuumTimeInt_;
}

const word& decoupledTimeData::couplingOption() const
{
    return couplingOption_;
}

const timeInterval& decoupledTimeData::molecularTimeInterval() const
{
    return molecularTimeInt_;
}

const bool& decoupledTimeData::continuumToMolecularBCs() const
{
    return continuumToMolecularBCs_;
}

const bool& decoupledTimeData::molecularToContinuumBCs() const
{
    return molecularToContinuumBCs_;    
}

// const bool& decoupledTimeData::molecularSwitch() const
// {
//     return molecularSwitch_;    
// }


//- place this function right at the end of the hybrid solver while loop
void decoupledTimeData::solversToRun()
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

//     if(couplingOption_ == "synchronous")
//     {
//         if
//         (
//             (!molecularTimeInt_.endTime()) && 
//             (continuumTimeInt_.endTime())
//         )
//         {
//             continuumSolution_ = false;
//             molecularSolution_ = true;
//         }
// 
//         // - just in case 
//         if
//         (
//             (molecularTimeInt_.endTime()) && 
//             (!continuumTimeInt_.endTime())
//         )
//         {
//             continuumSolution_ = true;
//             molecularSolution_ = false;
//         }
// 
//         if
//         (
//             (molecularTimeInt_.endTime()) && 
//             (continuumTimeInt_.endTime())
//         )
//         {
//             if(runSequential_) // continuum "block" first
//             {
//                 continuumSolution_ = true;
//                 molecularSolution_ = false;
//             }
//             else // continuum and molecular MD time-steps proceed
//             {
//                 continuumSolution_ = true;
//                 molecularSolution_ = true;
//             }
// 
//             continuumTimeInt_.endTime() = false;
//             molecularTimeInt_.endTime() = false;
//         }
//     }
}


void decoupledTimeData::boundaryConditions()
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

//     if(couplingOption_ == "synchronous")
//     {
//         //- boundary conditions tested here
//         if
//         (
//             (molecularTimeInt_.endTime()) &&
//             (continuumTimeInt_.endTime()) 
//         )
//         {
//             molecularToContinuumBCs_ = true;
//             continuumToMolecularBCs_ = true;
//         }
//         else
//         {
//             molecularToContinuumBCs_ = false;
//             continuumToMolecularBCs_ = false;
//         }
//     }
// 
//     if(couplingOption_ == "decoupled") // same as a sequential scheme
//     {
//         //- boundary conditions tested here
//         if(molecularTimeInt_.endTime())
//         {
//             molecularToContinuumBCs_ = true;
//         }
//         else
//         {
//             molecularToContinuumBCs_ = false;
//         }
//     
//         if(continuumTimeInt_.endTime())
//         {
//             continuumToMolecularBCs_ = true;
//         }
//         else
//         {
//             continuumToMolecularBCs_ = false;
//         }
//     }
}

// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

//- Prefix increment
decoupledTimeData& decoupledTimeData::operator++()
{
    if(molecularSolution_)
    {
        molecularTimeInt_++;
    }

    if(continuumSolution_)
    {
        continuumTimeInt_++;
    }

    boundaryConditions();

    return *this;
}


//- Postfix increment
decoupledTimeData& decoupledTimeData::operator++(int)
{
    return operator++();
}



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
