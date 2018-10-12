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
    molConfig

Description

SourceFiles
    molConfigI.H
    molConfig.C
    molConfigIO.C

\*---------------------------------------------------------------------------*/

#ifndef pointInfo_H
#define pointInfo_H

#include "scalar.H"
#include "vector.H"
#include "point.H"
#include "label.H"
#include "labelList.H"
#include "labelField.H"
#include "pTraits.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


class polyMesh;

/*---------------------------------------------------------------------------*\
                           Class pointInfo Declaration
\*---------------------------------------------------------------------------*/

class pointInfo
{

    // Private data

        point originPoint_;

        point currentPoint_;

        label currentPointLabel_;

        scalar Rsqr_;

        labelList patchList_;

        label startingProcNo_;

        bool zonePoint_;

        bool acceptedPoint_;

        bool visitedPoint_;

        bool boundaryPoint_;

        bool procBoundaryPoint_;

        bool cyclicBoundaryPoint_;

        bool interfacePoint_;

        bool switchedProc_;




    // Private Member Functions

        inline void setOriginPoint(const point& p);

        inline void setCurrentPoint(const point& p);

        inline void setCurrentPointLabel(const label& p);

        inline void setRsqr(const scalar& R);


//         //- Disallow default bitwise assignment
//         void operator=(const pointInfo&);

//         //- Disallow default bitwise copy construct
//         pointInfo(const pointInfo&);


public:

    // Static data members

        //- initial point far away.
        static point greatPoint;

    // Constructors

        //- Construct null
        inline pointInfo();


        //- Construct from currentPoint and currentPointLabel.
//         inline pointInfo
//         (
//            const point& currentPoint,
//            const label& currentPointLabel
//         );


        //- Construct from components: interfacePoint_, currentPoint_, Rsqr_
//         inline pointInfo
//         (
//             const point& originPoint,
//             const point& currentPoint,
//             const label& currentPointLabel,
//             const scalar& Rsqr
//         );


//         - Construct as copy
        inline pointInfo(const pointInfo&);


        //- Construct and return a clone
//         autoPtr<pointInfo> clone() const
//         {
//             return autoPtr<pointInfo>(new pointInfo(*this));
//         }

    // Destructor

        ~pointInfo();

    // Member Functions

        inline void setPatchList(const labelList& );

        inline void setPointInfo
        (
            const point& currentPoint,
            const label& currentPointLabel
        );


        // update the information of a point and test it
        inline void updatePointInfo
        (
            const point& originPoint,
            const scalar& Rsqr
        );

        // update the information of an interfacePoint
        inline void updateStartingPointInfo
        (
            const point& originPoint,
            const scalar& Rsqr
        );

        inline void updateStartingProcNo(const label& procN);

        // test within range
        inline void testPoint();


        // set functions

        inline void pointInZone();

        inline void pointOutOfZone();

        inline void acceptPoint();

        inline void resetAcceptedPoint();

        inline void visitPoint();

        inline void unVisitPoint();

        inline void pointOnBoundary();

        inline void pointOnProcBoundary();

        inline void pointOnCyclicBoundary();

        inline void pointOnInterface();

        inline void resetInterfacePoint();

        inline void switchProc();

        inline void unSwitchProc();


        // Access

            inline const point& originPoint() const;

            inline const point& currentPoint() const;

            inline const label& currentPointLabel() const;

            inline const scalar& Rsqr() const;

            inline const label& procNo() const;

            inline const labelList& patchList() const;

            inline const bool& isZonePoint() const;

            inline const bool& isAcceptedPoint() const;

//             inline const bool& isRejectedPoint() const;

            inline const bool& isVisitedPoint() const;

            inline const bool& isBoundaryPoint() const;

            inline const bool& isProcBoundaryPoint() const;

            inline const bool& isCyclicBoundaryPoint() const;

            inline const bool& isInterfacePoint() const;

            inline const bool& isSwitchedProc() const;


    // this is probably needed in order for the algorithm not to oscillate continously
    // between one processor and the next. I think we can do without it, but we must be 
    // careful of the implementation.



        // Check


        // Edit


        // Write

    // Member Operators

        //Note: Used to determine whether to call update.
//         inline bool operator==(const pointInfo&) const;
//
//         inline bool operator!=(const pointInfo&) const;

        inline friend bool operator==
        (
            const pointInfo& a,
            const pointInfo& b
        );

        inline friend bool operator!=
        (
            const pointInfo& a,
            const pointInfo& b
        );

    // IOstream Operators

        friend Ostream& operator<<(Ostream&, const pointInfo&);
        friend Istream& operator>>(Istream&, pointInfo&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "pointInfoI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
