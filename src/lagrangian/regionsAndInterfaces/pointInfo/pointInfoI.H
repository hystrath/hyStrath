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

Description

\*---------------------------------------------------------------------------*/

#include "polyMesh.H"

namespace Foam
{


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Null constructor
inline pointInfo::pointInfo()
:
    originPoint_(greatPoint),
    currentPoint_(greatPoint),
    currentPointLabel_(-1),
    Rsqr_(GREAT),
    patchList_(),
    startingProcNo_(-1),
    zonePoint_(false),
    acceptedPoint_(false),
//     rejectedPoint_(false),
    visitedPoint_(false),
    boundaryPoint_(false),
    procBoundaryPoint_(false),
    cyclicBoundaryPoint_(false),
    interfacePoint_(false),
    switchedProc_(false)
{}


// Construct from a currentPoint

// inline pointInfo::pointInfo
// (
//     const point& currentPoint,
//     const label& currentPointLabel
// )
// :
//     originPoint_(greatPoint),
//     currentPoint_(currentPoint),
//     currentPointLabel_(currentPointLabel),
//     Rsqr_(GREAT),
//     zonePoint_(false),
// //     visitedPoint_(true),
//     acceptedPoint_(false),
//     boundaryPoint_(false),
//     interfacePoint_(false)
// //     switchedDomain_(false)
// {}


// Construct from an interfacePoint, currentPoint, currentPointLabel and  Rsqr
// inline pointInfo::pointInfo
// (
//     const point& originPoint,
//     const point& currentPoint,
//     const label& currentPointLabel,
//     const scalar& Rsqr
// )
// :
//     originPoint_(originPoint),
//     currentPoint_(currentPoint),
//     currentPointLabel_(currentPointLabel),
//     Rsqr_(Rsqr),
//     zonePoint_(false),
// //     visitedPoint_(true),
//     acceptedPoint_(false),
//     boundaryPoint_(false),
//     interfacePoint_(false)
// //     switchedDomain_(false)
// {}


// Construct as copy
inline pointInfo::pointInfo(const pointInfo& wpt)
:
    originPoint_(wpt.originPoint()),
    currentPoint_(wpt.currentPoint()),
    currentPointLabel_(wpt.currentPointLabel()),
    Rsqr_(wpt.Rsqr()),
    patchList_(wpt.patchList()),
    startingProcNo_(wpt.procNo()),
    zonePoint_(wpt.isZonePoint()),
    acceptedPoint_(wpt.isAcceptedPoint()),
    visitedPoint_(wpt.isVisitedPoint()),
    boundaryPoint_(wpt.isBoundaryPoint()),
    procBoundaryPoint_(wpt.isProcBoundaryPoint()),
    cyclicBoundaryPoint_(wpt.isCyclicBoundaryPoint()),
    interfacePoint_(wpt.isInterfacePoint()),
    switchedProc_(wpt.isSwitchedProc())
{}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


inline void pointInfo::setOriginPoint(const point& p)
{
    originPoint_ = p;
}

inline void pointInfo::setCurrentPoint(const point& p)
{
    currentPoint_ = p;
}

inline void pointInfo::setCurrentPointLabel(const label& p)
{
    currentPointLabel_ = p;
}

inline void pointInfo::setRsqr(const scalar& R)
{
    Rsqr_ = R;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline void pointInfo::setPointInfo
(
    const point& currentPoint,
    const label& currentPointLabel
)
{
    setCurrentPoint(currentPoint);

    setCurrentPointLabel(currentPointLabel);
}


inline void pointInfo::setPatchList(const labelList& patchList)
{
    patchList_ = patchList;
}


inline void pointInfo::updateStartingPointInfo
(
    const point& originPoint,
    const scalar& Rsqr
)
{
    setOriginPoint(originPoint);

    setRsqr(Rsqr);

    acceptPoint();

    visitPoint();
}

inline void pointInfo::updateStartingProcNo(const label& procN)
{
   startingProcNo_ = procN;
}


// update the information of a point and test it
inline void pointInfo::updatePointInfo
(
    const point& originPoint,
    const scalar& Rsqr
)
{
    setOriginPoint(originPoint);

    setRsqr(Rsqr);


}

// test Point
inline void pointInfo::testPoint()
{
    scalar deltaRsqr = magSqr(originPoint_ - currentPoint_);

    if(deltaRsqr <= Rsqr_)
    {
        acceptPoint();
    }
}

// - Set

inline void pointInfo::pointInZone()
{
    zonePoint_ = true;
}

inline void pointInfo::pointOutOfZone()
{
    zonePoint_ = false;
}

inline void pointInfo::acceptPoint()
{
    acceptedPoint_ = true;
}

inline void pointInfo::resetAcceptedPoint()
{
    acceptedPoint_ = false;
}

inline void pointInfo::visitPoint()
{
    visitedPoint_ = true;
}

inline void pointInfo::unVisitPoint()
{
    visitedPoint_ = false;
}

inline void pointInfo::pointOnBoundary()
{
    boundaryPoint_ = true;
}

inline void pointInfo::pointOnProcBoundary()
{
    procBoundaryPoint_ = true;
}

inline void pointInfo::pointOnCyclicBoundary()
{
    cyclicBoundaryPoint_ = true;
}

inline void pointInfo::pointOnInterface()
{
    interfacePoint_ = true;
}

inline void pointInfo::resetInterfacePoint()
{
    interfacePoint_ = false;
}

inline void pointInfo::switchProc() 
{
    switchedProc_ = true;
}

inline void pointInfo::unSwitchProc() 
{
    switchedProc_ = false;
}

// - Access

inline const point& pointInfo::originPoint() const
{
    return originPoint_;
}

inline const point& pointInfo::currentPoint() const
{
    return currentPoint_;
}

inline const label& pointInfo::currentPointLabel() const
{
    return currentPointLabel_;
}

inline const scalar& pointInfo::Rsqr() const
{
    return Rsqr_;
}

inline const label& pointInfo::procNo() const
{
    return startingProcNo_;
}

inline const labelList& pointInfo::patchList() const
{
    return patchList_;
}

inline const bool& pointInfo::isZonePoint() const
{
    return zonePoint_;
}

inline const bool& pointInfo::isAcceptedPoint() const
{
    return acceptedPoint_;
}

inline const bool& pointInfo::isVisitedPoint() const
{
    return visitedPoint_;
}

inline const bool& pointInfo::isBoundaryPoint() const
{
    return boundaryPoint_;
}

inline const bool& pointInfo::isProcBoundaryPoint() const
{
    return procBoundaryPoint_;
}

inline const bool& pointInfo::isCyclicBoundaryPoint() const
{
    return cyclicBoundaryPoint_;
}

inline const bool& pointInfo::isInterfacePoint() const
{
    return interfacePoint_;
}

inline const bool& pointInfo::isSwitchedProc() const
{
    return switchedProc_;
}





// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //


inline bool operator==
(
    const pointInfo& a,
    const pointInfo& b
)
{
    return
    (
        a.originPoint() == b.originPoint()
     && a.currentPointLabel() == b.currentPointLabel()
     && a.currentPoint() == b.currentPoint()
     && a.Rsqr() == b.Rsqr()
     && a.isZonePoint()  == b.isZonePoint()
    );
}


inline bool operator!=
(
    const pointInfo& a,
    const pointInfo& b
)
{
    return !(a == b);
}


} // End namespace Foam

// ************************************************************************* //
