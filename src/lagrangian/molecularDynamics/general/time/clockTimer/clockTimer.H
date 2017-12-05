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
    clockTimer

Description

SourceFiles
    clockTimer.C

\*---------------------------------------------------------------------------*/

#ifndef clockTimer_H
#define clockTimer_H

#include "Time.H"
#include <stdio.h>
#include <ctime>
#include <unistd.h>
#include <sys/time.h>
#include "polyMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class clockTimer Declaration
\*---------------------------------------------------------------------------*/

class clockTimer
{
    // Private data

        Time& time_;

        word fieldName_;

//         timespec lastTime_;
        timeval lastTime_;
    
        scalar instTimeIndex_;
        scalar timeIndex_;
    
        // total duration time
        scalar totalDuration_;

        // duration in one MD time-step
        scalar duration_;
    
        // time over one write interval
        scalar instantDuration_;

        bool writeOut_;

        fileName timePath_;

public:

    // Constructors

        //- constructor
        clockTimer(Time&, word fieldName, bool write);

    // Destructor

        ~clockTimer();


    // Member functions

        void startClock();
        void stopClock();
        scalar averageTimeWriteInterval();
        scalar averageTime();
        void write();
        
        const scalar& instantDuration() const;
        const scalar& totalDuration() const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
