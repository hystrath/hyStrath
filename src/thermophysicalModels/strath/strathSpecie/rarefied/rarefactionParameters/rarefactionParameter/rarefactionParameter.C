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

#include "rarefactionParameter.H"
#include "dimensionedConstants.H"
#include "constants.H"

#include <string>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

    defineTypeNameAndDebug(rarefactionParameter, 0);
    defineRunTimeSelectionTable(rarefactionParameter, fvMesh);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::rarefactionParameter::rarefactionParameter
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

    computeRarefaction_(subDict("rarefiedParameters").lookupOrDefault<bool>("computeFieldAndBoundaries", true)),
    computeMfpBoundaries_(subDict("rarefiedParameters").lookupOrDefault<bool>("computeMfpBoundaries", true)),
    oldMfpDefinition_(subDict("rarefiedParameters").lookupOrDefault<bool>("oldMfpDefinition", false)),
    mfpModelName_(subDict("rarefiedParameters").lookup("mfpModel")),

    mfpMix_
    (
        IOobject
        (
            "mfp",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("mfp", dimLength, 0.0)
    ),

    Knov_
    (
        IOobject
        (
            "Kn_ov",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("Kn_ov", dimless, 0.0)
    ),

    characteristicLength_(readScalar(subDict("rarefiedParameters").lookup("characteristicLength"))),

    KnGLL_
    (
        IOobject
        (
            "KnGLL",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("KnGLL", dimless, 0.0)
    ),

    writeMfpSpecies_(subDict("rarefiedParameters").lookup("writeMfpSpecies", false)),
    writeMfpMixture_(subDict("rarefiedParameters").lookup("writeMfpMixture", false)),
    writeKnGLL_(subDict("rarefiedParameters").lookupOrDefault<bool>("writeKnGLL", false)),
    writeKnGLLComponents_(subDict("rarefiedParameters").lookupOrDefault<bool>("writeKnGLL_components", false)),
    writeKnOv_(subDict("rarefiedParameters").lookupOrDefault<bool>("writeKn_overall", false))

{
    if(thermo.hyLight())
    {
        computeRarefaction_ = false;
        writeMfpSpecies_ = false;
        writeMfpMixture_ = false;
        writeKnGLL_ = false;
        writeKnGLLComponents_ = false;
        writeKnOv_ = false;
    }

    const word dictThermoPhy
    (
        fileName(thermo.lookup("foamChemistryThermoFile")).name()
    );

    // Construct the mean free path model
    mfpModel_.set
    (
        new mfpModel
        (
            IOdictionary::name(),
            dictThermoPhy,
            species(),
            thermo.p(),
            thermo.T()
        )
    );

    mfp_.setSize(species().size());
    forAll(mfp_, speciei)
    {
        mfp_.set
        (
            speciei,
            new volScalarField
            (
                IOobject
                (
                    "mfp_" + species()[speciei],
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh_,
                dimensionedScalar("mfp_" + species()[speciei], dimLength, 0.0)
            )
        );
    }

    wordList KnsGLLNames(3, word::null);
    KnsGLLNames[0] = "rho"; KnsGLLNames[1] = "T"; KnsGLLNames[2] = "U";

    KnsGLL_.setSize(3);
    forAll(KnsGLL_, i)
    {
        KnsGLL_.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    "KnGLL_" + KnsGLLNames[i],
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh_,
                dimensionedScalar("KnGLL_" + KnsGLLNames[i], dimless, 0.0)
            )
        );
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
