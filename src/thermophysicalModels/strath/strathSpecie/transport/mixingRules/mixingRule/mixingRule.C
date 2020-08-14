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

\*---------------------------------------------------------------------------*/

#include "mixingRule.H"
#include "dimensionedConstants.H"
#include "constants.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

    defineTypeNameAndDebug(mixingRule, 0);
    defineRunTimeSelectionTable(mixingRule, fvMesh);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::mixingRule::mixingRule
(
    rho2ReactionThermo& thermo,
    const compressibleTurbulenceModel& turbulence
)
:
    IOdictionary
    (
        thermo.transportDictionary()
    ),

    mesh_(thermo.T().mesh()),
    thermo_(thermo),
    turbulence_(turbulence),

    spmu_(species().size()),
    spkappatr_(species().size()),
    spkappave_(species().size()),
    spalphatr_(species().size()),
    spalphave_(species().size()),

    writeMuSpecies_(subDict("transportModels").lookupOrDefault<bool>("writeViscositySpecies", false)),
    writeMuMixture_(subDict("transportModels").lookupOrDefault<bool>("writeViscosityMixture", false)),
    writeKappaSpecies_(subDict("transportModels").lookupOrDefault<bool>("writeThermalConducSpecies", false)),
    writeKappaMixture_(subDict("transportModels").lookupOrDefault<bool>("writeThermalConducMixture", false))

{
    forAll(spmu_, speciei)
    {
        spmu_.set
        (
            speciei,
            new volScalarField
            (
                IOobject
                (
                    "mu_" + species()[speciei],
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh_,
                dimensionedScalar("mu_" + species()[speciei], dimMass/dimLength/dimTime, 0.0)
            )
        );

        spkappatr_.set
        (
            speciei,
            new volScalarField
            (
                IOobject
                (
                    "kappatr_" + species()[speciei],
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh_,
                dimensionedScalar("kappatr_" + species()[speciei], dimensionSet(1,1,-3,-1,0,0,0), 0.0)
            )
        );

        spkappave_.set
        (
            speciei,
            new volScalarField
            (
                IOobject
                (
                    "kappave_" + species()[speciei],
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh_,
                dimensionedScalar("kappave_" + species()[speciei], dimensionSet(1,1,-3,-1,0,0,0), 0.0)
            )
        );

        spalphatr_.set
        (
            speciei,
            new volScalarField
            (
                IOobject
                (
                    "alphatr_" + species()[speciei],
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh_,
                dimensionedScalar("alphatr_" + species()[speciei], dimMass/dimLength/dimTime, 0.0)
            )
        );

        spalphave_.set
        (
            speciei,
            new volScalarField
            (
                IOobject
                (
                    "alphave_" + species()[speciei],
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh_,
                dimensionedScalar("alphave_" + species()[speciei], dimMass/dimLength/dimTime, 0.0)
            )
        );
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
