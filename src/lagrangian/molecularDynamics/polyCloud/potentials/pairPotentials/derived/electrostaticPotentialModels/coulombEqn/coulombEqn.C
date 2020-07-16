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

#include "coulombEqn.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(coulombEqn, 0);
addToRunTimeSelectionTable(pairPotentialModel, coulombEqn, dictionary);

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
coulombEqn::coulombEqn
(
    const polyMesh& mesh,
    polyMoleculeCloud& molCloud,
    const reducedUnits& redUnits,
    const word& name,
    const dictionary& dict
)
:
    pairPotentialModel(mesh, molCloud, redUnits, name, dict),
    constant_(1.0/(4.0 * constant::mathematical::pi * 8.854187817e-12))
{

    if(redUnits.runReducedUnits())
    {
        constant_ = (1.0/(4.0 * constant::mathematical::pi * redUnits.epsilonPermittivity()));
    }
    else
    {
        constant_ = 1.0/(4.0*constant::mathematical::pi*8.854187817e-12);
    }

    useTables_ = false;
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

coulombEqn::~coulombEqn()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

scalar coulombEqn::unscaledEnergy(const scalar r) const
{
    return constant_/r;
}

scalar coulombEqn::force(const scalar r) const
{
    scalar force = constant_*( 1/(r*r) );

    return force;
}

scalar coulombEqn::energy(const scalar r) const
{
    scalar energy = constant_*(1/r);

    return energy;
}

// void coulombEqn::interaction (const scalar r, scalar& force, scalar& energy)
// {
//     scalar oneOnR = (1/r);
//
//     force = constant_*(oneOnR*oneOnR);
//
//     energy = constant_*oneOnR;
// }

const dictionary& coulombEqn::dict() const
{
    return pairPotentialProperties_;
}

void  coulombEqn::write(const fileName& pathName)
{

}

} // End namespace Foam

// ************************************************************************* //
