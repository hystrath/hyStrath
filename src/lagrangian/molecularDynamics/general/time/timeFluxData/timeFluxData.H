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
    timeFluxData

Description
    Used as a time class framework for the controllers.


SourceFiles
    timeFluxData.C

\*---------------------------------------------------------------------------*/

#ifndef timeFluxData_H
#define timeFluxData_H

#include "Time.H"
#include "scalarField.H"
#include "timeInterval.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class timeFluxData Declaration
\*---------------------------------------------------------------------------*/

class timeFluxData
{
    // Private data

        Time& time_;

        scalar writeInterval_;
        label writeIntSteps_;


        timeInterval mdTime_;
        timeInterval averagingTime_;
        timeInterval controlTime_;
        timeInterval writeTime_;

        dimensionedScalar nAvTimeSteps_;

        // this is the ratio of number of averaging steps to the number of control steps (or their time-intervals)
        scalar nControlSteps_;

        label totalNAvSteps_;
        label totalNContSteps_;
        label totalNWrSteps_;

        label controlTimeIndex_;
        label averagingTimeIndex_;
        label writeTimeIndex_;
//         scalarField controlTimes_;
//         scalarField averagingTimes_;

//         void checkConsistency();

public:

    // Constructors

        //- constructor from Time
        timeFluxData
        (
            Time& t
        );

        //- Construct from Time and timeDict
        timeFluxData
        (
            Time& t,
            const dictionary& timeDict
        );


    // Destructor

        ~timeFluxData();


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

        const bool& averagingTime() const;
        const bool& controlTime() const;

        //- return the number of ...
        const label& nControls() const;
        const label& nAverages() const;

        const dimensionedScalar& nAvTimeSteps() const;
        const scalar& nControlSteps() const;

        const label& totalNAvSteps() const;
        const label& totalNContSteps() const;
        const label& totalNWrSteps() const;

        const label& controlTimeIndex() const;
        const label& averagingTimeIndex() const;
        const label& writeTimeIndex() const;

        const timeInterval& mdTimeInterval() const;
        const timeInterval& averageTimeInterval() const;
        const timeInterval& controlTimeInterval() const;
        const timeInterval& writeTimeInterval() const;

        timeInterval& controlTimeInterval();
        timeInterval& averageTimeInterval();

//         const scalarField& controlTimes() const;
//         const scalarField& averagingTimes() const;

        scalarField writeTimes();
        scalarField controlTimes();
        scalarField averagingTimes();

    // Member operators

        //- Prefix increment
        timeFluxData& operator++();

        //- Postfix increment
        timeFluxData& operator++(int);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
