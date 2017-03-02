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

#include "rhoChemistryHTC.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::hTCModels::rhoChemistryHTC::rhoChemistryHTC
(
    const word& modelType,
    const fvMesh& mesh
)
:
    rhoHTCModel(modelType, mesh),
    chemistryPtr_(rhoChemistryModel::New(mesh))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::hTCModels::rhoChemistryHTC::~rhoChemistryHTC()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::rhoReactionThermo&
Foam::hTCModels::rhoChemistryHTC::thermo()
{
    return chemistryPtr_->thermo();
}


const Foam::rhoReactionThermo&
Foam::hTCModels::rhoChemistryHTC::thermo() const
{
    return chemistryPtr_->thermo();
}


Foam::tmp<Foam::volScalarField>
Foam::hTCModels::rhoChemistryHTC::rho() const
{
    return chemistryPtr_->thermo().rho();
}


// ************************************************************************* //
