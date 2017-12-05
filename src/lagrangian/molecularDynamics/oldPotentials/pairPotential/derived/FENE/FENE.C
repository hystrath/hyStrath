/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2008-2009 OpenCFD Ltd.
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

\*---------------------------------------------------------------------------*/

#include "FENE.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace pairPotentials
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(FENE, 0);

addToRunTimeSelectionTable
(
    pairPotential,
    FENE,
    dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

FENE::FENE
(
    const word& name,
    const reducedUnits& rU,
    const dictionary& pairPotentialProperties
)
:
    pairPotential(name, rU, pairPotentialProperties),
    FENECoeffs_(pairPotentialProperties.subDict(typeName + "Coeffs")),
    k_(readScalar(FENECoeffs_.lookup("k"))),
    rO_(readScalar(FENECoeffs_.lookup("rO")))
{

    if(rU.runReducedUnits())
    {
        rO_ /= rU.refLength();
        k_ /= rU.refEnergy()/(rU.refLength()*rU.refLength());
    }

    Info << "FENE properties. rO = " << rO_ << ", k = " << k_ << endl; 

    setLookupTables(rU);
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

scalar FENE::unscaledEnergy(const scalar r) const
{

    scalar ir2 = (r/rO_)*(r/rO_);

    return -0.5*k_*rO_*rO_*Foam::log(1.0-ir2);
}


bool FENE::read
(
    const dictionary& pairPotentialProperties,
    const reducedUnits& rU
)
{
    pairPotential::read(pairPotentialProperties, rU);

    FENECoeffs_ = pairPotentialProperties.subDict(typeName + "Coeffs");

    FENECoeffs_.lookup("k") >> k_;
    FENECoeffs_.lookup("rO") >> rO_;

    if(rU.runReducedUnits())
    {
        rO_ /= rU.refLength();
        k_ /= rU.refEnergy()/(rU.refLength()*rU.refLength());
    }

    return true;
}


const dictionary& FENE::dict() const
{
    return FENECoeffs_;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace pairPotentials
} // End namespace Foam

// ************************************************************************* //
