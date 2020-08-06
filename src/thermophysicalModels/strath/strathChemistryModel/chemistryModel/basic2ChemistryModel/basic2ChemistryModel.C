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

#include "basic2ChemistryModel.H"
#include "fvMesh.H"
#include "Time.H"

/* * * * * * * * * * * * * * * private static data * * * * * * * * * * * * * */

namespace Foam
{
    defineTypeNameAndDebug(basic2ChemistryModel, 0);
}

// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

void Foam::basic2ChemistryModel::correct()
{
    // do nothing
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::basic2ChemistryModel::basic2ChemistryModel(const fvMesh& mesh)
:
    IOdictionary
    (
        IOobject
        (
            "chemistryProperties",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    ),
    mesh_(mesh),
    chemistry_(lookup("chemistry")),
    deltaTChemIni_(readScalar(lookup("initialChemicalTimeStep"))),
    deltaTChem_
    (
        IOobject
        (
            "deltaTChem",
            mesh.time().constant(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("deltaTChem0", dimTime, deltaTChemIni_)
    ),

    modifiedTemperature_(lookupOrDefault<Switch>("modifiedTemperature", false))

{
    if (modifiedTemperature_)
    {
        modTCoeffs_.setSize(2);
        modTCoeffs_.set
        (
            0,
            new scalar
            (
                subDict("modifiedTemperatureCoeffs")
                    .lookupOrDefault<scalar>("Tmin", 800.0)
            )
        );
        modTCoeffs_.set
        (
            1,
            new scalar
            (
                subDict("modifiedTemperatureCoeffs")
                    .lookupOrDefault<scalar>("epsilon", 80.0)
            )
        );
    }


    if (isDict("chemistryVibrationCoupling"))
    {
        CVModel_ = subDict("chemistryVibrationCoupling")
            .lookupOrDefault<word>("model", "none");
    }
    else
    {
        CVModel_ = "none";
    }

    if (CVModel_ == "ParkTTv")
    {
        exponentPark_ =
            subDict("chemistryVibrationCoupling").subDict(CVModel_ + "Coeffs")
                .lookupOrDefault<scalar>("exponentTtr", 0.7);

        ScvModel_ =
            subDict("chemistryVibrationCoupling").subDict(CVModel_ + "Coeffs")
                .lookupOrDefault<word>("sourceTermModel", "preferential");

        if (ScvModel_ != "preferential" and ScvModel_ != "nonPreferential")
        {
            FatalErrorIn
            (
                "Foam::basic2ChemistryModel::basic2ChemistryModel"
                "(const fvMesh& mesh)"
            )   << "Park's TTV sourceTerm model can either be 'preferential' "
                << "or 'nonPreferential'"
                << exit(FatalError);
        }
        else if (ScvModel_ == "preferential")
        {
            preferentialModel_ =
                subDict("chemistryVibrationCoupling")
                    .subDict(CVModel_ + "Coeffs").subDict("preferentialModel")
                    .lookupOrDefault<word>("factorType", word::null);

            if
            (
                preferentialModel_ != "constant" 
             && preferentialModel_ != "lineFitted"
            )
            {
                FatalErrorIn
                (
                    "Foam::basic2ChemistryModel::basic2ChemistryModel"
                    "(const fvMesh& mesh)"
                )   << "Park's TTV preferential model can either have "
                    << "'constant' or 'lineFitted' coefficients"
                    << exit(FatalError);
            }
        }
    }
    else if (CVModel_ != "CVDV" and CVModel_ != "none")
    {
        FatalErrorIn
        (
            "Foam::basic2ChemistryModel::basic2ChemistryModel(const fvMesh& "
            "mesh)"
        )
            << "This chemistry-vibration coupling model does not exist. "
            << "Valid models are 'ParkTTv', 'CVDV' and 'none'"
            << exit(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::basic2ChemistryModel::~basic2ChemistryModel()
{}


// ************************************************************************* //
