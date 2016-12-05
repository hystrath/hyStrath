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

#include "mixed2VELMixEnergyFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "multi2Thermo.H" // NEW VINCENT


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::mixed2VELMixEnergyFvPatchScalarField::
mixed2VELMixEnergyFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(p, iF)
{
    valueFraction() = 0.0;
    refValue() = 0.0;
    refGrad() = 0.0;
}


Foam::mixed2VELMixEnergyFvPatchScalarField::
mixed2VELMixEnergyFvPatchScalarField
(
    const mixed2VELMixEnergyFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchScalarField(ptf, p, iF, mapper)
{}


Foam::mixed2VELMixEnergyFvPatchScalarField::
mixed2VELMixEnergyFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchScalarField(p, iF, dict)
{}


Foam::mixed2VELMixEnergyFvPatchScalarField::
mixed2VELMixEnergyFvPatchScalarField
(
    const mixed2VELMixEnergyFvPatchScalarField& tppsf
)
:
    mixedFvPatchScalarField(tppsf)
{}


Foam::mixed2VELMixEnergyFvPatchScalarField::
mixed2VELMixEnergyFvPatchScalarField
(
    const mixed2VELMixEnergyFvPatchScalarField& tppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(tppsf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::mixed2VELMixEnergyFvPatchScalarField::updateCoeffs()
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
        const_cast<fvPatchScalarField&>(multiThermo.Tv().boundaryField()[patchi])
    );

    Tvw.evaluate();

    valueFraction() = Tvw.valueFraction();
    refValue() = multiThermo.hevel(pw, Tvw.refValue(), patchi);
    refGrad() =
        multiThermo.Cv_v(pw, Tvw, patchi)*Tvw.refGrad()
      + patch().deltaCoeffs()*
        (
            multiThermo.hevel(pw, Tvw, patchi)
          - multiThermo.hevel(pw, Tvw, patch().faceCells())
        );

    mixedFvPatchScalarField::updateCoeffs();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        mixed2VELMixEnergyFvPatchScalarField
    );
}

// ************************************************************************* //
