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
    
    modifiedTemperature_(lookupOrDefault<Switch>("modifiedTemperature", false)) // NEW VINCENT 16/02/2017
    
{
    if(modifiedTemperature_)
    {
        modTCoeffs_.setSize(2);
        modTCoeffs_.set(0, new scalar(readScalar(subDict("modifiedTemperatureCoeffs").lookup("Tmin"))));
        modTCoeffs_.set(1, new scalar(readScalar(subDict("modifiedTemperatureCoeffs").lookup("epsilon"))));
    }
    
    
    if(isDict("chemistryVibrationCoupling")) // NEW VINCENT 13/08/2016
    {
        CVModel_ = subDict("chemistryVibrationCoupling").lookupOrDefault<word>("model", word::null);
    }
    else
    {
        CVModel_ = "none";
    }
    
    if (CVModel_ == "ParkTTv")
    {
        exponentPark_ = readScalar(subDict("chemistryVibrationCoupling").subDict(CVModel_ + "Coeffs").lookup("exponentTtr"));
        
        ScvModel_ = subDict("chemistryVibrationCoupling").subDict(CVModel_ + "Coeffs").lookupOrDefault<word>("sourceTermModel", "preferential");
        
        if (ScvModel_ != "preferential" and ScvModel_ != "nonPreferential")
        {
            FatalErrorIn("Foam::basic2ChemistryModel::basic2ChemistryModel(const fvMesh& mesh)")
                << "ScvModel_ model can either be 'preferential' or 'nonPreferential'."
                << exit(FatalError);
        }
        
        if (ScvModel_ == "preferential")
        {
            preferentialModel_ = subDict("chemistryVibrationCoupling").subDict(CVModel_ + "Coeffs").subDict("preferentialModel").lookupOrDefault<word>("factorType", word::null);
            
            if (preferentialModel_ != "constant" and preferentialModel_ != "lineFitted")
            {
                FatalErrorIn("Foam::basic2ChemistryModel::basic2ChemistryModel(const fvMesh& mesh)")
                    << "ScvModel_ preferential model can either be defined as 'constant' or 'lineFitted'."
                    << exit(FatalError);
            }    
        } 
    }
    else if (CVModel_ != "CVDV" and CVModel_ != "none")
    {
        FatalErrorIn("Foam::basic2ChemistryModel::basic2ChemistryModel(const fvMesh& mesh)")
            << "This chemistryVibrationCoupling model does not exist. Valid models are 'ParkTTv' and 'CVDV'."
            << exit(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::basic2ChemistryModel::~basic2ChemistryModel()
{}


// ************************************************************************* //
