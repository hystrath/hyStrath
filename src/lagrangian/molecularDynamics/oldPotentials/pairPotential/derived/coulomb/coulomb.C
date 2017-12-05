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

#include "coulomb.H"
#include "mathematicalConstants.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace pairPotentials
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(coulomb, 0);

addToRunTimeSelectionTable
(
    pairPotential,
    coulomb,
    dictionary
);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

coulomb::coulomb
(
    const word& name,
    const reducedUnits& rU,
    const dictionary& pairPotentialProperties
)
:
    pairPotential(name, rU, pairPotentialProperties),
    oneOverFourPiEps0_(1.0/(4.0 * constant::mathematical::pi * 8.854187817e-12))
{
    if(rU.runReducedUnits())
    {
        oneOverFourPiEps0_ = (1.0/(4.0 * constant::mathematical::pi * rU.epsilonPermittivity()));
    }
    else
    {
    	oneOverFourPiEps0_ = 1.0/(4.0*constant::mathematical::pi*8.854187817e-12);
    }

    setLookupTables(rU);
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

scalar coulomb::unscaledEnergy(const scalar r) const
{
    return oneOverFourPiEps0_/r;
}


bool coulomb::read
(
    const dictionary& pairPotentialProperties,
    const reducedUnits& rU
)
{
    pairPotential::read(pairPotentialProperties, rU);

    return true;
}

const dictionary& coulomb::dict() const
{
    return pairPotentialProperties_;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace pairPotentials
} // End namespace Foam

// ************************************************************************* //
