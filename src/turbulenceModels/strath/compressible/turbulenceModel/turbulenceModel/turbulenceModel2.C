/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
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

#include "turbulenceModel2.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "fvcGrad.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(turbulenceModel2, 0);
defineRunTimeSelectionTable(turbulenceModel2, turbulenceModel2);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

turbulenceModel2::turbulenceModel2
(
    const volScalarField& rho,
    const volVectorField& U,
    const surfaceScalarField& phi,
    const rho2ReactionThermo& thermophysicalModel,
    const word& turbulenceModel2Name
)
:
    regIOobject
    (
        IOobject
        (
            turbulenceModel2Name,
            U.time().constant(),
            U.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        )
    ),
    runTime_(U.time()),
    mesh_(U.mesh()),

    rho_(rho),
    U_(U),
    phi_(phi),
    thermophysicalModel_(thermophysicalModel),

    y_(mesh_)
{}


// * * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * //

autoPtr<turbulenceModel2> turbulenceModel2::New
(
    const volScalarField& rho,
    const volVectorField& U,
    const surfaceScalarField& phi,
    const rho2ReactionThermo& thermophysicalModel,
    const word& turbulenceModel2Name
)
{
    // get model name, but do not register the dictionary
    // otherwise it is registered in the database twice
    const word modelType
    (
        IOdictionary
        (
            IOobject
            (
                "turbulenceProperties",
                U.time().constant(),
                U.db(),
                IOobject::MUST_READ_IF_MODIFIED,
                IOobject::NO_WRITE,
                false
            )
        ).lookup("simulationType")
    );

    Info<< "Selecting turbulence model type " << modelType << endl;

    turbulenceModel2ConstructorTable::iterator cstrIter =
        turbulenceModel2ConstructorTablePtr_->find(modelType);

    if (cstrIter == turbulenceModel2ConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "turbulenceModel2::New(const volScalarField&, "
            "const volVectorField&, const surfaceScalarField&, "
            "rho2ReactionThermo&, const word&)"
        )   << "Unknown turbulenceModel2 type "
            << modelType << nl << nl
            << "Valid turbulenceModel2 types:" << endl
            << turbulenceModel2ConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<turbulenceModel2>
    (
        cstrIter()(rho, U, phi, thermophysicalModel, turbulenceModel2Name)
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

tmp<volScalarField> turbulenceModel2::rhoEpsilonEff() const
{
    tmp<volTensorField> tgradU = fvc::grad(U_);
    return mu()*(tgradU() && dev(twoSymm(tgradU()))) + rho_*epsilon();
}


void turbulenceModel2::correct()
{
    if (mesh_.changing())
    {
        y_.correct();
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace compressible
} // End namespace Foam

// ************************************************************************* //
