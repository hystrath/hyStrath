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

#include "mixed2VELEnergyFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "multi2Thermo.H" // NEW VINCENT

#include <string.H> // NEW VINCENT 18/04/2016

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::mixed2VELEnergyFvPatchScalarField::
mixed2VELEnergyFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(p, iF),
    thermo_(rho2ReactionThermo::lookup2ReactionThermo(*this))
{
    word fieldName = iF.name();
    specieName_ = fieldName.substr(fieldName.find("_") + 1);
    
    valueFraction() = 0.0;
    refValue() = 0.0;
    refGrad() = 0.0;
}


Foam::mixed2VELEnergyFvPatchScalarField::
mixed2VELEnergyFvPatchScalarField
(
    const mixed2VELEnergyFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchScalarField(ptf, p, iF, mapper),
    thermo_(rho2ReactionThermo::lookup2ReactionThermo(*this))
{
    word fieldName = iF.name();
    specieName_ = fieldName.substr(fieldName.find("_") + 1);
}


Foam::mixed2VELEnergyFvPatchScalarField::
mixed2VELEnergyFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchScalarField(p, iF, dict),
    thermo_(rho2ReactionThermo::lookup2ReactionThermo(*this))
{
    word fieldName = iF.name();
    specieName_ = fieldName.substr(fieldName.find("_") + 1);
}


Foam::mixed2VELEnergyFvPatchScalarField::
mixed2VELEnergyFvPatchScalarField
(
    const mixed2VELEnergyFvPatchScalarField& tppsf
)
:
    mixedFvPatchScalarField(tppsf),
    thermo_(rho2ReactionThermo::lookup2ReactionThermo(*this))
{}


Foam::mixed2VELEnergyFvPatchScalarField::
mixed2VELEnergyFvPatchScalarField
(
    const mixed2VELEnergyFvPatchScalarField& tppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(tppsf, iF),
    thermo_(rho2ReactionThermo::lookup2ReactionThermo(*this))
{
    word fieldName = iF.name();
    specieName_ = fieldName.substr(fieldName.find("_") + 1);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::mixed2VELEnergyFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const multi2Thermo& multiThermo = multi2Thermo::lookup2Thermo(*this);
    const label patchi = patch().index();

    const scalarField& pw = multiThermo.p().boundaryField()[patchi];
    mixedFvPatchScalarField& Tvw = refCast<mixedFvPatchScalarField>
    (
        const_cast<fvPatchScalarField&>(thermo_.composition().Tv(specieName_).boundaryField()[patchi])
    );

    Tvw.evaluate();

    valueFraction() = Tvw.valueFraction();
    refValue() = thermo_.composition().hevel(specieName_, pw, Tvw.refValue(), patchi);
    refGrad() =
        thermo_.composition().Cv_vel(specieName_, pw, Tvw, patchi)*Tvw.refGrad()
      + patch().deltaCoeffs()*
        (
            thermo_.composition().hevel(specieName_, pw, Tvw, patchi)
          - thermo_.composition().hevel(specieName_, pw, Tvw, patch().faceCells())
        );

    mixedFvPatchScalarField::updateCoeffs();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        mixed2VELEnergyFvPatchScalarField
    );
}

// ************************************************************************* //
