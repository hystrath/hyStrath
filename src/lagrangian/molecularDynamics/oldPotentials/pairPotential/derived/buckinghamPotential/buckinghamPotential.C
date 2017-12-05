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

#include "buckinghamPotential.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace pairPotentials
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(buckinghamPotential, 0);

addToRunTimeSelectionTable
(
    pairPotential,
    buckinghamPotential,
    dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

buckinghamPotential::buckinghamPotential
(
    const word& name,
    const reducedUnits& rU,
    const dictionary& pairPotentialProperties
)
:
    pairPotential(name, rU, pairPotentialProperties),
    buckinghamPotentialCoeffs_(pairPotentialProperties.subDict(typeName + "Coeffs")),
    rho_(readScalar(buckinghamPotentialCoeffs_.lookup("rho"))),
    alpha_(readScalar(buckinghamPotentialCoeffs_.lookup("alpha"))),
    C_(readScalar(buckinghamPotentialCoeffs_.lookup("C")))
{

    if(rU.runReducedUnits())
    {
        rho_ /= rU.refLength();
        alpha_ /= rU.refEnergy();
        
        C_ /= rU.refEnergy();

        scalar l = rU.refLength();
        scalar lSix = l*l*l*l*l*l;

        C_ /= lSix;
    }


    setLookupTables(rU);
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

scalar buckinghamPotential::unscaledEnergy(const scalar r) const
{
    scalar ir2 = (1.0/r)*(1.0/r);

    scalar ir6 = ir2*ir2*ir2;

    return alpha_*exp(-r/rho_) - C_*ir6;
}


bool buckinghamPotential::read
(
    const dictionary& pairPotentialProperties,
    const reducedUnits& rU
)
{
    pairPotential::read(pairPotentialProperties, rU);

    buckinghamPotentialCoeffs_ = pairPotentialProperties.subDict(typeName + "Coeffs");

    buckinghamPotentialCoeffs_.lookup("rho") >> rho_;
    buckinghamPotentialCoeffs_.lookup("alpha") >> alpha_;
    buckinghamPotentialCoeffs_.lookup("C") >> C_;

    if(rU.runReducedUnits())
    {
        rho_ /= rU.refLength();
        alpha_ /= rU.refEnergy();
        
        C_ /= rU.refEnergy();

        scalar l = rU.refLength();
        scalar lSix = l*l*l*l*l*l;

        C_ /= lSix;
    }

    return true;
}


const dictionary& buckinghamPotential::dict() const
{
    return buckinghamPotentialCoeffs_;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace pairPotentials
} // End namespace Foam

// ************************************************************************* //
