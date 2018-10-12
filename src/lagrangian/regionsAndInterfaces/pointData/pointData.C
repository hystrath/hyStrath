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

#include "pointData.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

Foam::point Foam::pointData::greatPoint(GREAT, GREAT, GREAT);



// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os, const pointData& wDist)
{
    if (os.format() == IOstream::ASCII)
    {
        return os
            << wDist.origin() << token::SPACE << wDist.distSqr()
            << token::SPACE << wDist.s() << token::SPACE << wDist.v();
    }
    else
    {
        return os
            << wDist.origin() << wDist.distSqr() << wDist.s() << wDist.v();
    }
}

Foam::Istream& Foam::operator>>(Istream& is, pointData& wDist)
{
    return is >> wDist.origin_ >> wDist.distSqr_ >> wDist.s_ >> wDist.v_;
}

// ************************************************************************* //
