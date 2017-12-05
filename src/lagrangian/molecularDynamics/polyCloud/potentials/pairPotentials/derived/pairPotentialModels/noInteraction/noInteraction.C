/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2007 OpenCFD Ltd.
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

Description

\*---------------------------------------------------------------------------*/

#include "noInteraction.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(noInteraction, 0);
addToRunTimeSelectionTable(pairPotentialModel, noInteraction, dictionary);

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
noInteraction::noInteraction
(
    const polyMesh& mesh,
    polyMoleculeCloud& molCloud,
    const reducedUnits& redUnits,
    const word& name, 
    const dictionary& dict
)
:
    pairPotentialModel(mesh, molCloud, redUnits, name, dict) 
{
//     exclusions_ = true;
    useTables_ = false;
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

noInteraction::~noInteraction()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

scalar noInteraction::unscaledEnergy(const scalar r) const
{
    return 0.0;
}


scalar noInteraction::force(const scalar r) const
{
//     return forceLookUpFromTable(r);
    return 0.0;
}
    
scalar noInteraction::energy(const scalar r) const
{
//     return energyLookUpFromTable(r);
    return 0.0;
}


// bool noInteraction::read
// (
//     const dictionary& pairPotentialProperties,
//     const reducedUnits& rU
// )
// {
//     pairPotentialModel::read(pairPotentialProperties, rU);
// 
//     noInteractionCoeffs_ = pairPotentialProperties.subDict(typeName + "Coeffs");
// 
//     noInteractionCoeffs_.lookup("sigma") >> sigma_;
//     noInteractionCoeffs_.lookup("epsilon") >> epsilon_;
// 
//     if(rU.runReducedUnits())
//     {
//         sigma_ /= rU.refLength();
//         epsilon_ /= rU.refEnergy();
//     }
// 
//     return true;
// }

const dictionary& noInteraction::dict() const
{
    return propsDict_;
}
void noInteraction::write(const fileName& pathName)
{
    
}

} // End namespace Foam

// ************************************************************************* //
