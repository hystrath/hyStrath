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

Class
    timeData

Description
    Used as a time class framework for the controllers.


SourceFiles
    timeData.C

\*---------------------------------------------------------------------------*/

#ifndef timeData_H
#define timeData_H

#include "Time.H"
#include "scalarField.H"
#include "timeInterval.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class timeData Declaration
\*---------------------------------------------------------------------------*/

class timeData
{
    // Private data

        Time& time_;

        dictionary timeDict_;

        scalar writeInterval_;
        label writeIntSteps_;

        bool resetFieldsAtOutput_;
        label resetIndex_;

        timeInterval mdTime_;
        timeInterval samplingTime_;
        timeInterval averagingTime_;
//         timeInterval calcPropTime_;
        timeInterval controlTime_;

        dimensionedScalar nAvTimeSteps_;

        // this is the ratio of number of averaging steps to the number of control steps (or their time-intervals)
        scalar nControlSteps_;

        label totalNSampSteps_;
        label totalNAvSteps_;
        label totalNContSteps_;
//         label totalNCalcSteps_;

        label controlTimeIndex_;
//         label calcTimeIndex_;
        label averagingTimeIndex_;

        bool decoupledFromWriteInterval_;

//         scalarField controlTimes_;
//         scalarField calcTimes_;
//         scalarField averagingTimes_;

//         void checkConsistency();

public:

    // Constructors

        //- constructor from Time
        timeData
        (
            Time& t
        );

        //- Construct from Time and timeDict
        timeData
        (
            Time& t,
            const dictionary& timeDict
        );


    // Destructor

        ~timeData();


    // Member functions
        void setInitialData();

        void checkAndModifyTimeProperties();

        void setTimeData(const dictionary& timeDict);


        //-return constant access to time
        const Time& time() const;
        Time& time();

        //- write interval
        const scalar& writeInterval() const;
        const label& writeIntervalSteps() const;


        const bool& samplingTime() const;
        const bool& averagingTime() const;
        const bool& controlTime() const;
//         const bool& calcPropTime() const;

        //- return the number of ...
        const label& nSamples() const;
        const label& nControls() const;
        const label& nAverages() const;
//         const label& nCalcProp() const;

        const dimensionedScalar& nAvTimeSteps() const;
        const scalar& nControlSteps() const;

        // for resetting
        scalar nAveragingTimeSteps();


        const label& totalNSampSteps() const;
        const label& totalNAvSteps() const;
        const label& totalNContSteps() const;
//         const label& totalNCalcSteps() const;

        const label& controlTimeIndex() const;
//         const label& calcTimeIndex() const;
        const label& averagingTimeIndex() const;

        const timeInterval& mdTimeInterval() const;
        const timeInterval& sampleTimeInterval() const;
        const timeInterval& averageTimeInterval() const;
//         const timeInterval& calcPropTimeInterval() const;
        const timeInterval& controlTimeInterval() const;
        
        timeInterval& controlTimeInterval();
//         timeInterval& calcPropTimeInterval();
        timeInterval& averageTimeInterval();

//         const scalarField& controlTimes() const;
//         const scalarField& calcTimes() const;
//         const scalarField& averagingTimes() const;


        const bool& resetFieldsAtOutput() const;
        bool& resetFieldsAtOutput();

        scalarField controlTimes();
        scalarField averagingTimes();




    // Member operators

        //- Prefix increment
        timeData& operator++();

        //- Postfix increment
        timeData& operator++(int);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
