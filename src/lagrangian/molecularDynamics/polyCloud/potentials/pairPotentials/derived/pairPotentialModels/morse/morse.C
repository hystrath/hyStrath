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

Description

\*---------------------------------------------------------------------------*/

#include "morse.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(morse, 0);
addToRunTimeSelectionTable(pairPotentialModel, morse, dictionary);

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
morse::morse
(
    const polyMesh& mesh,
    polyMoleculeCloud& molCloud,
    const reducedUnits& redUnits,
    const word& name,
    const dictionary& dict
)
:
    pairPotentialModel(mesh, molCloud, redUnits, name, dict),
    propsDict_(dict.subDict(typeName + "Coeffs")),
    Kcr_(readScalar(propsDict_.lookup("Kcr"))),
    gamma_(readScalar(propsDict_.lookup("gamma"))),
    rC_(readScalar(propsDict_.lookup("rC")))
{
    if(redUnits.runReducedUnits())
    {
        Kcr_ /= redUnits.refEnergy();
        gamma_ *= redUnits.refLength();
        rC_ /= redUnits.refLength();
    }

    setLookupTables();
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

morse::~morse()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

scalar morse::unscaledEnergy(const scalar r) const
{
    scalar exponent = -gamma_*(r-rC_);
    scalar exp = Foam::exp(exponent);

    return Kcr_*(exp-1.0)*(exp-1.0);
}

scalar morse::force(const scalar r) const
{
    return forceLookUpFromTable(r);
}

scalar morse::energy(const scalar r) const
{
    return energyLookUpFromTable(r);
}


// bool morse::read
// (
//     const dictionary& pairPotentialProperties,
//     const reducedUnits& rU
// )
// {
//     pairPotentialModel::read(pairPotentialProperties, rU);
//
//     morseCoeffs_ = pairPotentialProperties.subDict(typeName + "Coeffs");
//
//     morseCoeffs_.lookup("sigma") >> sigma_;
//     morseCoeffs_.lookup("epsilon") >> epsilon_;
//
//     if(rU.runReducedUnits())
//     {
//         sigma_ /= rU.refLength();
//         epsilon_ /= rU.refEnergy();
//     }
//
//     return true;
// }
void morse::write(const fileName& pathName)
{

}

const dictionary& morse::dict() const
{
    return propsDict_;
}


} // End namespace Foam

// ************************************************************************* //
