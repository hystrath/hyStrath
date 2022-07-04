/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016-2021 hyStrath
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
    trackingNumber

Description

\*----------------------------------------------------------------------------*/

#include "trackingNumber.H"

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::trackingNumber::resetTrackingNumbers()
{
    label trackingNumber = trackingIndex_;

    //- parallel-processing
    if(Pstream::parRun())
    {
        trackingNumber = trackingIndex_ + Pstream::myProcNo() + (Pstream::nProcs()-1)*trackingIndex_;

        //- sending
        for (label p = 0; p < Pstream::nProcs(); p++)
        {
            if(p != Pstream::myProcNo())
            {
                const label proc = p;
                {
                    OPstream toNeighbour(Pstream::commsTypes::blocking, proc);
                    toNeighbour << trackingNumber;
                }
            }
        }

        //- receiving
        for (label p = 0; p < Pstream::nProcs(); p++)
        {
            if(p != Pstream::myProcNo())
            {
                label trackingNumberProc;

                const label proc = p;
                {
                    IPstream fromNeighbour(Pstream::commsTypes::blocking, proc);
                    fromNeighbour >> trackingNumberProc;
                }

                if(trackingNumberProc > trackingNumber)
                {
                    trackingNumber = trackingNumberProc;
                }
            }
        }
    }

//     label maxLimit = GREAT;

    if(trackingNumber >= maxLimit_)
    {
        resetTracking_ = true;
        trackingIndex_ = 0;

        Info << "WARNING : resetting tracking numbers." << endl;
    }
    else
    {
        resetTracking_ = false;
    }
}


Foam::label Foam::trackingNumber::getTrackingNumber()
{
    label trackingNumber = trackingIndex_;

    if (Pstream::parRun())
    {
        //note: for parallel processing the available set of tracking numbers gets
        // divided according to the number of processors being used.

        trackingNumber = trackingIndex_ + Pstream::myProcNo() + (Pstream::nProcs()-1)*trackingIndex_;
    }

    trackingIndex_++;

    return trackingNumber;
}

Foam::label Foam::trackingNumber::getMaxTrackingNumber()
{
    label trackingNumber = trackingIndex_;

    if (Pstream::parRun())
    {
        trackingNumber = trackingIndex_ + Pstream::myProcNo() + (Pstream::nProcs()-1)*trackingIndex_;

        //- sending
        for (label p = 0; p < Pstream::nProcs(); p++)
        {
            if(p != Pstream::myProcNo())
            {
                const label proc = p;
                {
                    OPstream toNeighbour(Pstream::commsTypes::blocking, proc);
                    toNeighbour << trackingNumber;
                }
            }
        }

        //- receiving
        for (label p = 0; p < Pstream::nProcs(); p++)
        {
            if(p != Pstream::myProcNo())
            {
                label trackingNumberProc;

                const label proc = p;
                {
                    IPstream fromNeighbour(Pstream::commsTypes::blocking, proc);
                    fromNeighbour >> trackingNumberProc;
                }

                if(trackingNumberProc > trackingNumber)
                {
                    trackingNumber = trackingNumberProc;
                }
            }
        }
    }

    return trackingNumber;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Null constructor
trackingNumber::trackingNumber()
:
    trackingIndex_(0),
    resetTracking_(false),
    maxLimit_(labelMax)
{
//     Info << "max: " << labelMax << endl;

}



// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

trackingNumber::~trackingNumber()
{}


Foam::label& Foam::trackingNumber::trackingIndex()
{
    return trackingIndex_;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //



// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //



// * * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
