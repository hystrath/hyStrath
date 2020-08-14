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

#include "relaxationTimeModelVV.H"
#include "dimensionedConstants.H"
#include "constants.H"

#include <string>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

    defineTypeNameAndDebug(relaxationTimeModelVV, 0);
    defineRunTimeSelectionTable(relaxationTimeModelVV, fvMesh);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::relaxationTimeModelVV::relaxationTimeModelVV
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
    const word dict2T(IOdictionary::name()), dictThermoPhy
    (
        fileName(thermo.lookup("foamChemistryThermoFile")).name()
    );

    // Construct the relaxation time model
    tauVVijModel_.set
    (
        new VVModel
        (
            dict2T,
            dictThermoPhy,
            solvedVibEqSpecies(), // NEW VINCENT 06/08/2016
            thermo.p(),
            thermo.T(),
            thermo.composition().Tv(),
            thermo.composition().nD()
        )
    );

    QVV_.setSize(solvedVibEqSpecies().size());


    forAll(solvedVibEqSpecies(), speciei)
    {
        QVV_.set
        (
            speciei,
            new volScalarField
            (
                IOobject
                (
                    "QVV_" + solvedVibEqSpecies()[speciei],
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh_,
                dimensionedScalar("QVV", dimensionSet(1,-1,-3,0,0), 0.0)
            )
        );
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

/*Foam::tmp<Foam::volScalarField>
Foam::relaxationTimeModelVV::VVRelaxationSource()
{
    tmp<volScalarField> tQVV
    (
        new volScalarField
        (
            IOobject
            (
                "VVRelaxationSource",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar("QVV", dimensionSet(1, -1, -3, 0, 0), 0.0)
        )
    );

    return tQVV;
}*/


bool Foam::relaxationTimeModelVV::read()
{
    return regIOobject::read();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
