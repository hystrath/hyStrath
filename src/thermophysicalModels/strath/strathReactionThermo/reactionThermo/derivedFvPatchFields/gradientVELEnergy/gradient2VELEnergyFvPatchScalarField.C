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

#include "gradient2VELEnergyFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"

#include "multi2Thermo.H" // NEW VINCENT 20/02/2016
#include <string.H> // NEW VINCENT 20/02/2016

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::gradient2VELEnergyFvPatchScalarField::
gradient2VELEnergyFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(p, iF),
    thermo_(rho2ReactionThermo::lookup2ReactionThermo(*this))
{
    word fieldName = iF.name();
    specieName_ = fieldName.substr(fieldName.find("_") + 1);
} // Only this constructor is used at run-time


Foam::gradient2VELEnergyFvPatchScalarField::
gradient2VELEnergyFvPatchScalarField
(
    const gradient2VELEnergyFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedGradientFvPatchScalarField(ptf, p, iF, mapper),
    thermo_(rho2ReactionThermo::lookup2ReactionThermo(*this))
{
    word fieldName = iF.name();
    specieName_ = fieldName.substr(fieldName.find("_") + 1);
}


Foam::gradient2VELEnergyFvPatchScalarField::
gradient2VELEnergyFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedGradientFvPatchScalarField(p, iF, dict),
    thermo_(rho2ReactionThermo::lookup2ReactionThermo(*this))
{
    word fieldName = iF.name();
    specieName_ = fieldName.substr(fieldName.find("_") + 1);
}


Foam::gradient2VELEnergyFvPatchScalarField::
gradient2VELEnergyFvPatchScalarField
(
    const gradient2VELEnergyFvPatchScalarField& tppsf
)
:
    fixedGradientFvPatchScalarField(tppsf),
    thermo_(rho2ReactionThermo::lookup2ReactionThermo(*this))
{}


Foam::gradient2VELEnergyFvPatchScalarField::
gradient2VELEnergyFvPatchScalarField
(
    const gradient2VELEnergyFvPatchScalarField& tppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(tppsf, iF),
    thermo_(rho2ReactionThermo::lookup2ReactionThermo(*this))
{
    word fieldName = iF.name();
    specieName_ = fieldName.substr(fieldName.find("_") + 1);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::gradient2VELEnergyFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }
    
    //Info << "gradient2VELEnergy is used for patch called " << patch().name() << ", species " << specieName_ << endl;

    const multi2Thermo& multiThermo = multi2Thermo::lookup2Thermo(*this);
    const label patchi = patch().index();

    const scalarField& pw = multiThermo.p().boundaryField()[patchi];
    
    fvPatchScalarField& spTvw =
        const_cast<fvPatchScalarField&>(thermo_.composition().Tv(specieName_).boundaryField()[patchi]);
    spTvw.evaluate();
    
    gradient() = thermo_.composition().Cv_vel(specieName_, pw, spTvw, patchi)*spTvw.snGrad()
      + patch().deltaCoeffs()*
        (
            thermo_.composition().hevel(specieName_, pw, spTvw, patchi)
          - thermo_.composition().hevel(specieName_, pw, spTvw, patch().faceCells())
        );
        
    fixedGradientFvPatchScalarField::updateCoeffs();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        gradient2VELEnergyFvPatchScalarField
    );
}

// ************************************************************************* //
