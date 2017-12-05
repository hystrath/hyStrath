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

#include "coulomb.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(coulomb, 0);
addToRunTimeSelectionTable(pairPotentialModel, coulomb, dictionary);

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
coulomb::coulomb
(
    const polyMesh& mesh,
    polyMoleculeCloud& molCloud,
    const reducedUnits& redUnits,
    const word& name, 
    const dictionary& dict
)
:
    pairPotentialModel(mesh, molCloud, redUnits, name, dict),
    oneOverFourPiEps0_(1.0/(4.0 * constant::mathematical::pi * 8.854187817e-12))   
{
 
    if(redUnits.runReducedUnits())
    {
        oneOverFourPiEps0_ = (1.0/(4.0 * constant::mathematical::pi * redUnits.epsilonPermittivity()));
    }
    else
    {
        oneOverFourPiEps0_ = 1.0/(4.0*constant::mathematical::pi*8.854187817e-12);
    }

    setLookupTables();   
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

coulomb::~coulomb()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

scalar coulomb::unscaledEnergy(const scalar r) const
{
    return oneOverFourPiEps0_/r;
}



scalar coulomb::force(const scalar r) const
{
    return forceLookUpFromTable(r);
}
    
scalar coulomb::energy(const scalar r) const
{
    return energyLookUpFromTable(r);
}

// void coulombEqn::interaction (const scalar r, scalar& force, scalar& energy)
// {
//     force = forceLookUpFromTable(r);
//     
//     energy = energyLookUpFromTable(r);
// }

const dictionary& coulomb::dict() const
{
    return pairPotentialProperties_;
}

void  coulomb::write(const fileName& pathName)
{
    
}

} // End namespace Foam

// ************************************************************************* //
