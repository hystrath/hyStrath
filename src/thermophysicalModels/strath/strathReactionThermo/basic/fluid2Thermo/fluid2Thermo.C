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

#include "fluid2Thermo.H"

/* * * * * * * * * * * * * * * private static data * * * * * * * * * * * * * */

namespace Foam
{
    defineTypeNameAndDebug(fluid2Thermo, 0);
    defineRunTimeSelectionTable(fluid2Thermo, fvMesh);
}

const Foam::word Foam::fluid2Thermo::transportDictName("transportProperties");


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fluid2Thermo::fluid2Thermo(const fvMesh& mesh, const word& phaseName)
:
    basic2Thermo(mesh, phaseName),
    
    transportDictionary_
    (
        IOobject
        (
            phasePropertyName(transportDictName, phaseName),
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    )
{}


Foam::fluid2Thermo::fluid2Thermo
(
    const fvMesh& mesh,
    const dictionary& dict,
    const word& phaseName
)
:
    basic2Thermo(mesh, dict, phaseName),
    
    transportDictionary_
    (
        IOobject
        (
            phasePropertyName(transportDictName, phaseName),
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    )
{}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::fluid2Thermo> Foam::fluid2Thermo::New
(
    const fvMesh& mesh,
    const word& phaseName
)
{
    return basic2Thermo::New<fluid2Thermo>(mesh, phaseName);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fluid2Thermo::~fluid2Thermo()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::fluid2Thermo::nu() const
{
    return mu()/rho();
}


Foam::tmp<Foam::scalarField> Foam::fluid2Thermo::nu(const label patchi) const
{
    return mu(patchi)/rho(patchi);
}

// ************************************************************************* //
