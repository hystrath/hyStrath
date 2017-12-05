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

#include "morse.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace pairPotentials
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(morse, 0);

addToRunTimeSelectionTable
(
    pairPotential,
    morse,
    dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

morse::morse
(
    const word& name,
    const reducedUnits& rU,
    const dictionary& pairPotentialProperties
)
:
    pairPotential(name, rU, pairPotentialProperties),
    morseCoeffs_(pairPotentialProperties.subDict(typeName + "Coeffs")),
    Kcr_(readScalar(morseCoeffs_.lookup("Kcr"))),
    gamma_(readScalar(morseCoeffs_.lookup("gamma"))),
    rC_(readScalar(morseCoeffs_.lookup("rC")))    
{

    if(rU.runReducedUnits())
    {
        Kcr_ /= rU.refEnergy();
        gamma_ *= rU.refLength();        
        rC_ /= rU.refLength();
    }

    setLookupTables(rU);
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

scalar morse::unscaledEnergy(const scalar r) const
{
    scalar exponent = -gamma_*(r-rC_);
    scalar exp = Foam::exp(exponent);
    
    return Kcr_*(exp-1.0)*(exp-1.0);
}


bool morse::read
(
    const dictionary& pairPotentialProperties,
    const reducedUnits& rU
)
{
    pairPotential::read(pairPotentialProperties, rU);

    morseCoeffs_ = pairPotentialProperties.subDict(typeName + "Coeffs");

    morseCoeffs_.lookup("Kcr") >> Kcr_;
    morseCoeffs_.lookup("gamma") >> gamma_;
    morseCoeffs_.lookup("rC") >> rC_;    

    if(rU.runReducedUnits())
    {
        Kcr_ /= rU.refEnergy();
        gamma_ *= rU.refLength();        
        rC_ /= rU.refLength();
    }

    return true;
}


const dictionary& morse::dict() const
{
    return morseCoeffs_;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace pairPotentials
} // End namespace Foam

// ************************************************************************* //
