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

#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "multi2Thermo.H" // NEW VINCENT
#include "addToRunTimeSelectionTable.H"
#include "fixed2VELEnergyFvPatchScalarField.H"

#include <string.H> // NEW VINCENT 20/02/2016

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fixed2VELEnergyFvPatchScalarField::
fixed2VELEnergyFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF),
    thermo_(rho2ReactionThermo::lookup2ReactionThermo(*this))
{
    word fieldName = iF.name();
    specieName_ = fieldName.substr(fieldName.find("_") + 1);
} // Only this constructor is used at run-time


Foam::fixed2VELEnergyFvPatchScalarField::
fixed2VELEnergyFvPatchScalarField
(
    const fixed2VELEnergyFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper),
    thermo_(rho2ReactionThermo::lookup2ReactionThermo(*this))
{
    word fieldName = iF.name();
    specieName_ = fieldName.substr(fieldName.find("_") + 1);
}


Foam::fixed2VELEnergyFvPatchScalarField::
fixed2VELEnergyFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF, dict),
    thermo_(rho2ReactionThermo::lookup2ReactionThermo(*this))
{
    word fieldName = iF.name();
    specieName_ = fieldName.substr(fieldName.find("_") + 1);
}


Foam::fixed2VELEnergyFvPatchScalarField::
fixed2VELEnergyFvPatchScalarField
(
    const fixed2VELEnergyFvPatchScalarField& tppsf
)
:
    fixedValueFvPatchScalarField(tppsf),
    thermo_(rho2ReactionThermo::lookup2ReactionThermo(*this))
{}


Foam::fixed2VELEnergyFvPatchScalarField::
fixed2VELEnergyFvPatchScalarField
(
    const fixed2VELEnergyFvPatchScalarField& tppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(tppsf, iF),
    thermo_(rho2ReactionThermo::lookup2ReactionThermo(*this))
{
    word fieldName = iF.name();
    specieName_ = fieldName.substr(fieldName.find("_") + 1);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fixed2VELEnergyFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    //Info << "fixed2VELEnergy is used for patch called " << patch().name() << ", species " << specieName_ << endl;
    
    const multi2Thermo& multiThermo = multi2Thermo::lookup2Thermo(*this);
    
    const label patchi = patch().index();

    const scalarField& pw = multiThermo.p().boundaryField()[patchi];

    fvPatchScalarField& spTvw =
        const_cast<fvPatchScalarField&>(thermo_.composition().Tv(specieName_).boundaryField()[patchi]);
    spTvw.evaluate();
        
    operator==(thermo_.composition().hevel(specieName_, pw, spTvw, patchi));

    fixedValueFvPatchScalarField::updateCoeffs();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        fixed2VELEnergyFvPatchScalarField
    );
}

// ************************************************************************* //
