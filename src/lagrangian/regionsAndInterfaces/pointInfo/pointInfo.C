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
    pointInfo

Description

\*----------------------------------------------------------------------------*/

#include "pointInfo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

point pointInfo::greatPoint(GREAT, GREAT, GREAT);


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

Ostream& operator<<(Ostream& os, const pointInfo& p)
{
    if (os.format() == IOstream::ASCII)
    {
        return os
            << p.currentPointLabel() << p.originPoint() <<  p.Rsqr();
    }
    else
    {
        return os
            <<  p.currentPointLabel() << p.originPoint() <<  p.Rsqr();
    }
}

Istream& operator>>(Istream& is, pointInfo& p)
{
    return is >> p.currentPointLabel_ >> p.originPoint_ >>  p.Rsqr_;
}


//             inline const point& originPoint() const;
//
//             inline const point& currentPoint() const;
//
//             inline const label& currentPointLabel() const;
//
//             inline const scalar& Rsqr() const;
//
//
//             inline const bool& isZonePoint() const;
//
//             inline const bool& isAcceptedPoint() const;
//
// //             inline const bool& isRejectedPoint() const;
//
//             inline const bool& isVisitedPoint() const;
//
//             inline const bool& isBoundaryPoint() const;
//
//             inline const bool& isProcBoundaryPoint() const;
//
//             inline const bool& isCyclicBoundaryPoint() const;
//
//             inline const bool& isInterfacePoint() const;

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

pointInfo::~pointInfo()
{}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //




// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
