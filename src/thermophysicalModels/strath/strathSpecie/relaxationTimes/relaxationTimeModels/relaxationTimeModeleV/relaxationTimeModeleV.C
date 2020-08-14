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

#include "relaxationTimeModeleV.H"
#include "dimensionedConstants.H"
#include "constants.H"

#include <string>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

    defineTypeNameAndDebug(relaxationTimeModeleV, 0);
    defineRunTimeSelectionTable(relaxationTimeModeleV, fvMesh);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::relaxationTimeModeleV::relaxationTimeModeleV
(
    rho2ReactionThermo& thermo,
    const compressibleTurbulenceModel& turbulence
)
:
    IOdictionary
    (
        thermo.twoTemperatureDictionary()
    ),

    mesh_(thermo.T().mesh()),
    thermo_(thermo),
    turbulence_(turbulence)

{
    const word dictThermoPhy
    (
        fileName(thermo.lookup("foamChemistryThermoFile")).name()
    );

    // Construct the relaxation time model
    taueViModel_.set
    (
        new eVModel
        (
            IOdictionary::name(),
            dictThermoPhy,
            species(),
            thermo.composition().pP("e-"),
            thermo.composition().Tv("e-")
        )
    );

    QeV_.setSize(solvedVibEqSpecies().size());

    forAll(solvedVibEqSpecies(), speciei)
    {
        QeV_.set
        (
            speciei,
            new volScalarField
            (
                IOobject
                (
                    "QeV_" + species()[speciei],
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh_,
                dimensionedScalar("QeV", dimensionSet(1,-1,-3,0,0), 0.0)
            )
        );
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

/*Foam::tmp<Foam::volScalarField>
Foam::relaxationTimeModeleV::eVRelaxationSource()
{
    tmp<volScalarField> tQeV
    (
        new volScalarField
        (
            IOobject
            (
                "eVRelaxationSource",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar("QeV", dimensionSet(1, -1, -3, 0, 0), 0.0)
        )
    );

    return tQeV;
}*/


bool Foam::relaxationTimeModeleV::read()
{
    return regIOobject::read();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
