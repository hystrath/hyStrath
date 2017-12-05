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
    timeDataMeas

Description

SourceFiles
    timeDataMeas.C

\*---------------------------------------------------------------------------*/

#ifndef timeDataMeas_H
#define timeDataMeas_H

#include "Time.H"
#include "scalarField.H"
#include "timeInterval.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class timeDataMeas Declaration
\*---------------------------------------------------------------------------*/

class timeDataMeas
{
    // Private data

        Time& time_;

        dictionary timeDict_;

        word timeMeasOption_;

        scalar writeInterval_;
        label writeIntSteps_;

        bool resetFieldsAtOutput_;
        label resetIndex_;

        timeInterval mdTime_;
        timeInterval samplingTime_;
        timeInterval averagingTime_;
        timeInterval writeTime_;

//         timeInterval calcPropTime_;

        dimensionedScalar nAvTimeSteps_;

        label totalNSampSteps_;
        label totalNAvSteps_;
//         label totalNCalcSteps_;

        label averagingTimeIndex_;
        label samplingTimeIndex_;

//         scalarField averagingTimes_;
//         scalarField samplingTimes_;

//         void checkConsistency();

public:

    // Constructors

        //- constructor from Time
        timeDataMeas
        (
            Time& t
        );

        //- Construct from Time and timeDict
        timeDataMeas
        (
            Time& t,
            const dictionary& timeDict
        );


    // Destructor

        ~timeDataMeas();


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
        const bool& writeTime() const;
//         const bool& calcPropTime() const;

        //- return the number of ...
        const label& nSamples() const;
        const label& nAverages() const;
//         const label& nCalcProp() const;

        const dimensionedScalar& nAvTimeSteps() const;

        const bool& resetFieldsAtOutput() const;
        bool& resetFieldsAtOutput();

        // for resetting
        scalar nAveragingTimeSteps();



        const label& totalNSampSteps() const;
        const label& totalNAvSteps() const;
//         const label& totalNCalcSteps() const;

        const label& averagingTimeIndex() const;
        const label& samplingTimeIndex() const;
//         const label& calcTimeIndex() const;

        const timeInterval& mdTimeInterval() const;
        const timeInterval& sampleTimeInterval() const;
        const timeInterval& averageTimeInterval() const;
        const timeInterval& writeTimeInterval() const;


        const word& timeOption() const;

//         const timeInterval& calcPropTimeInterval() const;
       
//         timeInterval& calcPropTimeInterval();

//         const scalarField& averagingTimes() const;
//         const scalarField& samplingTimes() const;
        scalarField averagingTimes();
        scalarField samplingTimes();
        scalarField writeTimes();
        label writeSteps();



    // Member operators

        //- Prefix increment
        timeDataMeas& operator++();

        //- Postfix increment
        timeDataMeas& operator++(int);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
