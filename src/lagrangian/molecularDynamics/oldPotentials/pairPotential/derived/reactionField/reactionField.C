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

#include "reactionField.H"
#include "mathematicalConstants.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace pairPotentials
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(reactionField, 0);

addToRunTimeSelectionTable
(
    pairPotential,
    reactionField,
    dictionary
);



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

reactionField::reactionField
(
    const word& name,
    const reducedUnits& rU,
    const dictionary& pairPotentialProperties
)
:
    pairPotential(name, rU, pairPotentialProperties),
    reactionFieldCoeffs_
    (
        pairPotentialProperties.subDict(typeName + "Coeffs")
    ),
    dielectricConst_(readScalar(reactionFieldCoeffs_.lookup("dielectricConst"))), 
    oneOverFourPiEps0_(1.0/(4.0 * constant::mathematical::pi * 8.854187817e-12))
{
    if(rU.runReducedUnits())
    {
        oneOverFourPiEps0_ = (1.0/(4.0 * constant::mathematical::pi * rU.epsilonPermittivity()));
    }

    setLookupTables(rU);
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

scalar reactionField::unscaledEnergy(const scalar r) const
{
    scalar B = (2.0*(dielectricConst_ - 1.0))/((2.0*dielectricConst_) + 1.0);

    scalar rCut = pairPotential::rCut_;

    scalar C = (3.0*dielectricConst_)/((2.0*dielectricConst_ + 1.0)*rCut);

    return oneOverFourPiEps0_*(((1.0/r) + ((B*(r*r))/(2*(rCut*rCut*rCut)))) - C );
}


bool reactionField::read
(
    const dictionary& pairPotentialProperties,
    const reducedUnits& rU
)
{
    pairPotential::read(pairPotentialProperties, rU);

    reactionFieldCoeffs_ =
        pairPotentialProperties.subDict(typeName + "Coeffs");

    reactionFieldCoeffs_.lookup("dielectricConst") >> dielectricConst_;

    return true;
}

const dictionary& reactionField::dict() const
{
    return reactionFieldCoeffs_;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace pairPotentials
} // End namespace Foam

// ************************************************************************* //
