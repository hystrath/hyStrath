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

#include "relaxationTimeModel.H"
#include "dimensionedConstants.H"
#include "constants.H"
#include <string>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

    defineTypeNameAndDebug(relaxationTimeModel, 0);
    defineRunTimeSelectionTable(relaxationTimeModel, fvMesh);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::relaxationTimeModel::relaxationTimeModel
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

    // Construct the V-T relaxation time model
    tauVTijModel_.set
    (
        new VTModel
        (
            dict2T,
            dictThermoPhy,
            molecules(),
            species(),
            thermo.p(),
            thermo.T(),
            thermo.composition().Tv(),
            thermo.composition().nD()
        )
    );

    QVT_.setSize(molecules().size());
    //QVTmode_.setSize(molecules().size()); // ABORTIVE WORK

    forAll(molecules(), moli)
    {
        QVT_.set
        (
            moli,
            new volScalarField
            (
                IOobject
                (
                    "QVT_" + molecules()[moli],
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh_,
                dimensionedScalar("QVT", dimensionSet(1, -1, -3, 0, 0), 0.0)
            )
        );
    }

    /*forAll(QVTmode_, moli) // ABORTIVE WORK
    {
        QVTmode_.set
        (
            moli,
            new PtrList<volScalarField>
            (
                thermo.composition().noVibrationalTemp(moli)
            )
        );
    }

    forAll(QVTmode_, moli)
    {
        forAll(QVTmode_[moli], vibMode)
        {
            QVTmode_[moli].set
            (
                vibMode,
                new volScalarField
                (
                    IOobject
                    (
                        "QVT_" + molecules()[moli] + "." + word(vibMode+1),
                        mesh_.time().timeName(),
                        mesh_,
                        IOobject::NO_READ,
                        IOobject::NO_WRITE
                    ),
                    mesh_,
                    dimensionedScalar("QVT", dimensionSet(1, -1, -3, 0, 0), 0.0)
                )
            );
        }
    }*/
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::relaxationTimeModel::QVT()
{
    tmp<volScalarField> tQVT
    (
        new volScalarField
        (
            IOobject
            (
                "QVT",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh_,
            dimensionedScalar("QVT", dimensionSet(1, -1, -3, 0, 0), 0.0)
        )
    );
    
    scalarField& QVTCells = tQVT.ref().primitiveFieldRef();
    
    forAll(molecules(), moli)
    {
        forAll(QVTCells, celli)
        {
            QVTCells[celli] += QVT_[moli][celli];
        }
    }

    return tQVT;
}


bool Foam::relaxationTimeModel::read()
{
    return regIOobject::read();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
