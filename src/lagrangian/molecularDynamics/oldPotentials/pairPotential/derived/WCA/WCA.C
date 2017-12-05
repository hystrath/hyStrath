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

#include "WCA.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace pairPotentials
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(WCA, 0);

addToRunTimeSelectionTable
(
    pairPotential,
    WCA,
    dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

WCA::WCA
(
    const word& name,
    const reducedUnits& rU,
    const dictionary& pairPotentialProperties
)
:
    pairPotential(name, rU, pairPotentialProperties),
    WCACoeffs_(pairPotentialProperties.subDict(typeName + "Coeffs")),
    sigma_(readScalar(WCACoeffs_.lookup("sigma"))),
    epsilon_(readScalar(WCACoeffs_.lookup("epsilon")))
{

    if(rU.runReducedUnits())
    {
        sigma_ /= rU.refLength();
        epsilon_ /= rU.refEnergy();
    }

    setLookupTables(rU);
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

scalar WCA::unscaledEnergy(const scalar r) const
{
    // (rIJ/sigma)^-2
    scalar ir2 = (sigma_/r)*(sigma_/r);

    // (rIJ/sigma)^-6
    scalar ir6 = ir2*ir2*ir2;

    return 4.0 * epsilon_*(ir6*(ir6 - 1.0));
}


bool WCA::read
(
    const dictionary& pairPotentialProperties,
    const reducedUnits& rU
)
{
    pairPotential::read(pairPotentialProperties, rU);

    WCACoeffs_ = pairPotentialProperties.subDict(typeName + "Coeffs");

    WCACoeffs_.lookup("sigma") >> sigma_;
    WCACoeffs_.lookup("epsilon") >> epsilon_;

    if(rU.runReducedUnits())
    {
        sigma_ /= rU.refLength();
        epsilon_ /= rU.refEnergy();
    }

    return true;
}


const dictionary& WCA::dict() const
{
    return WCACoeffs_;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace pairPotentials
} // End namespace Foam

// ************************************************************************* //
