/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016-2020 hyStrath
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of hyStrath, a derivative work of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    timeInterval

Description

SourceFiles
    timeInterval.C

\*---------------------------------------------------------------------------*/

#ifndef timeInterval_H
#define timeInterval_H

#include "Time.H"
// #include "scalarField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class timeInterval Declaration
\*---------------------------------------------------------------------------*/

class timeInterval
{
    // Private data

        // - number of steps within timeInterval
        label nSteps_;

        //- time step
        scalar deltaT_;

        // - starts from zero or any number
//         label offsetTimeIndex_;

        // - continously changing time-index, linked to the runTime timer
        label timeIndex_;

        //- when timeIndex becomes equal to nSteps, this becomes true, otherwise
        // it is always false
        bool endTime_;

public:

    // Constructors

        //- Null constructor
        timeInterval();

        //- Construct from no of steps
        timeInterval
		(
            const label& nSteps
        );

        //- Construct from components
        timeInterval
        (
            const label& nSteps,
            const scalar& deltaT
        );

    // Destructor

        ~timeInterval();


    // Member functions

        const bool& endTime() const;
        bool& endTime();
        const label& nSteps() const;
        label& nSteps();
        const label& timeIndex() const;
        label& timeIndex();
        const scalar& deltaT() const;
        scalar& deltaT();
//         label& offsetTimeIndex();

    // Member operators

        //- Prefix increment
        timeInterval& operator++();

        //- Postfix increment
        timeInterval& operator++(int);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
