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
    decoupledTimeData

Description

SourceFiles
    decoupledTimeData.C

\*---------------------------------------------------------------------------*/

#ifndef decoupledTimeData_H
#define decoupledTimeData_H

#include "Time.H"
#include "scalarField.H"
#include "timeInterval.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class decoupledTimeData Declaration
\*---------------------------------------------------------------------------*/

class decoupledTimeData
{
    // Private data

        Time& continuumTime_;
        Time& molecularTime_;

        dictionary couplingTimeDict_;

        scalar contWriteInterval_;
        scalar molWriteInterval_;

        scalar contCouplingTimeInterval_;
        scalar molCouplingTimeInterval_;
//         label writeIntSteps_;

        timeInterval molecularTimeInt_;
        timeInterval continuumTimeInt_;

        bool continuumSolution_;
        bool molecularSolution_;

        bool molecularToContinuumBCs_;
        bool continuumToMolecularBCs_;


        // - sequential or synchronous
        word couplingOption_;
// 
//         bool decoupledScheme_;
//         bool molecularSwitch_;
//         
// 
//         bool runSequential_;


//         timeInterval calcPropTime_;

//         dimensionedScalar nAvTimeSteps_;

/*        label totalNSampSteps_;
        label totalNAvSteps_;*/
//         label totalNCalcSteps_;

//         label averagingTimeIndex_;
//         label samplingTimeIndex_;

//         scalarField averagingTimes_;
//         scalarField samplingTimes_;

//         void checkConsistency();

public:

    // Constructors

        //- constructor from Time
        decoupledTimeData
        (
            Time& continuumTime,
            Time& molecularTime
        );

        //- Construct from Time and timeDict
        decoupledTimeData
        (
            Time& continuumTime,
            Time& molecularTime,
            const dictionary& timeDict
        );


    // Destructor

        ~decoupledTimeData();


    // Member functions
//         void setInitialData();

//         void checkAndModifyTimeProperties();

//         void setTimeData(const dictionary& timeDict);

        //-return constant access to time
        const Time& continuumTime() const;
        Time& continuumTime();

        const Time& molecularTime() const;
        Time& molecularTime();

        //- write interval
        const scalar& contWriteInterval() const;
        const scalar& molWriteInterval() const;

        const timeInterval& continuumTimeInterval() const;
        const timeInterval& molecularTimeInterval() const;

        const bool& continuumSolution() const;
        bool& continuumSolution();

        const bool& molecularSolution() const;
        bool& molecularSolution();

        const word& couplingOption() const;

        const bool& continuumToMolecularBCs() const;
        const bool& molecularToContinuumBCs() const;


//         const bool& molecularSwitch() const;


        scalar couplingRatio();

        void solversToRun();

        void boundaryConditions();

    // Member operators

        //- Prefix increment
        decoupledTimeData& operator++();

        //- Postfix increment
        decoupledTimeData& operator++(int);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
