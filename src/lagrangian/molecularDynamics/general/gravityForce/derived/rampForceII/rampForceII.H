/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2005 OpenCFD Ltd.
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
   forceAtTimeZero : force at time t = 0
   forceAtRampEndTime : maximum ramp force
   rampEndTime : time at which maximum ramp force is reached
 
   optional
     relaxationTime : allows for an initial equilibration/relaxation time
     provided options to fixGradient (i.e. changes rampEndTime)
     or fix rampEndTime (and change gradient)

     acrossMultipleRuns : allows you to pick up where your simulation has finished.
     No need to change your other options.
 
Description

SourceFiles
    rampForceIII.H
    rampForceII.C
    rampForceIIIO.C

\*---------------------------------------------------------------------------*/

#ifndef rampForceII_H
#define rampForceII_H


#include "gravityForce.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class rampForceII Declaration
\*---------------------------------------------------------------------------*/

class rampForceII
:
    public gravityForce
{

private:

    // Private data

        dictionary propsDict_;

        vector force_;
        
        vector direction_;
        
        scalar initialForce_;
        
        scalar finalForce_;
        
        scalar endTime_;
        
        scalar forceGradient_;

        scalar relaxationTime_;
        
        scalar currentTime_;

        scalar deltaTMD_;
    
        scalar c_;
    
    

public:

    //- Runtime type information
    TypeName("rampForceII");

    // Static data members

    // Constructors

        //- construct from components
        rampForceII
        (
            Time& time,
            const dictionary& dict
        );

    // Selectors


    // Destructor

        ~rampForceII();

    // Member Functions

        //- output data

        // Access

        // Check

        // Edit

            void updateForce();

            vector force(const vector& position);
            
            vector force(const scalar& time);
            
        // Write


            void write
            (
                const fileName& fixedPathName,
                const fileName& timePath
            );


            void updateProperties
            (
                const dictionary& dict
            );


    // Member Operators

    // Friend Functions

    // Friend Operators

    // IOstream Operators

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
