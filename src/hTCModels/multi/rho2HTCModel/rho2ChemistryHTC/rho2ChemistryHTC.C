/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012 OpenFOAM Foundation
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

#include "rho2ChemistryHTC.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::hTC2Models::rho2ChemistryHTC::rho2ChemistryHTC
(
    const word& modelType,
    const fvMesh& mesh
)
:
    rho2HTCModel(modelType, mesh),
    chemistryPtr_(rho2ChemistryModel::New(mesh))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::hTC2Models::rho2ChemistryHTC::~rho2ChemistryHTC()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::rho2ReactionThermo&
Foam::hTC2Models::rho2ChemistryHTC::thermo()
{
    return chemistryPtr_->thermo();
}


const Foam::rho2ReactionThermo&
Foam::hTC2Models::rho2ChemistryHTC::thermo() const
{
    return chemistryPtr_->thermo();
}


Foam::tmp<Foam::volScalarField>
Foam::hTC2Models::rho2ChemistryHTC::rho() const
{
    return chemistryPtr_->thermo().rho();
}


// ************************************************************************* //
