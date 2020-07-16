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
    uniformForce

Description

\*----------------------------------------------------------------------------*/

#include "uniformForce.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(uniformForce, 0);

addToRunTimeSelectionTable(gravityForce, uniformForce, dictionary);



// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

//- Construct from components
uniformForce::uniformForce
(
    Time& time,
    const dictionary& dict
)
:
    gravityForce(time, dict),
    propsDict_(dict.subDict(typeName + "Properties")),
    force_(propsDict_.lookup("force"))
{
    timeVarying_ = true;
}



// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

uniformForce::~uniformForce()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

vector uniformForce::force(const vector& position)
{
    return vector::zero;
}

void uniformForce::updateForce()
{

}

vector uniformForce::force(const scalar& time)
{
    return force_;
}

void uniformForce::write
(
    const fileName& fixedPathName,
    const fileName& timePath
)
{

}

void uniformForce::updateProperties
(
    const dictionary& dict
)
{}

// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
