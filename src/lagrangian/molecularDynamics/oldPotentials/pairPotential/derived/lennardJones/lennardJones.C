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

#include "lennardJones.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace pairPotentials
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(lennardJones, 0);

addToRunTimeSelectionTable
(
    pairPotential,
    lennardJones,
    dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

lennardJones::lennardJones
(
    const word& name,
    const reducedUnits& rU,
    const dictionary& pairPotentialProperties
)
:
    pairPotential(name, rU, pairPotentialProperties),
    lennardJonesCoeffs_(pairPotentialProperties.subDict(typeName + "Coeffs")),
    sigma_(readScalar(lennardJonesCoeffs_.lookup("sigma"))),
    epsilon_(readScalar(lennardJonesCoeffs_.lookup("epsilon")))
{

    if(rU.runReducedUnits())
    {
        sigma_ /= rU.refLength();
        epsilon_ /= rU.refEnergy();
    }

    setLookupTables(rU);
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

scalar lennardJones::unscaledEnergy(const scalar r) const
{
    // (rIJ/sigma)^-2
    scalar ir2 = (sigma_/r)*(sigma_/r);

    // (rIJ/sigma)^-6
    scalar ir6 = ir2*ir2*ir2;

    return 4.0 * epsilon_*(ir6*(ir6 - 1.0));
}


bool lennardJones::read
(
    const dictionary& pairPotentialProperties,
    const reducedUnits& rU
)
{
    pairPotential::read(pairPotentialProperties, rU);

    lennardJonesCoeffs_ = pairPotentialProperties.subDict(typeName + "Coeffs");

    lennardJonesCoeffs_.lookup("sigma") >> sigma_;
    lennardJonesCoeffs_.lookup("epsilon") >> epsilon_;

    if(rU.runReducedUnits())
    {
        sigma_ /= rU.refLength();
        epsilon_ /= rU.refEnergy();
    }

    return true;
}


const dictionary& lennardJones::dict() const
{
    return lennardJonesCoeffs_;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace pairPotentials
} // End namespace Foam

// ************************************************************************* //
