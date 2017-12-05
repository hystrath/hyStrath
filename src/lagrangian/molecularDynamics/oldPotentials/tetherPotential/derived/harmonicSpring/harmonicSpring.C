/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

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

\*---------------------------------------------------------------------------*/

#include "harmonicSpring.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace tetherPotentials
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(harmonicSpring, 0);

addToRunTimeSelectionTable
(
    tetherPotential,
    harmonicSpring,
    dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

harmonicSpring::harmonicSpring
(
    const word& name,
    const reducedUnits& rU,
    const dictionary& tetherPotentialProperties
)
:
    tetherPotential(name, rU, tetherPotentialProperties),
    harmonicSpringCoeffs_
    (
        tetherPotentialProperties.subDict(typeName + "Coeffs")
    ),
    springConstant_(readScalar(harmonicSpringCoeffs_.lookup("springConstant")))
{
    if(rU.runReducedUnits())
    {
        springConstant_ /= (rU.refForce()/rU.refLength());
    }
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

scalar harmonicSpring::energy(const vector r) const
{
    return 0.5*springConstant_*magSqr(r);
}


vector harmonicSpring::force(const vector r) const
{
    return -springConstant_*r;
}


bool harmonicSpring::read
(
    const dictionary& tetherPotentialProperties,
    const reducedUnits& rU
)
{
    tetherPotential::read(tetherPotentialProperties, rU);

    harmonicSpringCoeffs_ =
        tetherPotentialProperties.subDict(typeName + "Coeffs");

    harmonicSpringCoeffs_.lookup("springConstant") >> springConstant_;

    if(rU.runReducedUnits())
    {
        springConstant_ /= (rU.refForce()/rU.refLength());
    }

    return true;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace tetherPotentials
} // End namespace Foam

// ************************************************************************* //
