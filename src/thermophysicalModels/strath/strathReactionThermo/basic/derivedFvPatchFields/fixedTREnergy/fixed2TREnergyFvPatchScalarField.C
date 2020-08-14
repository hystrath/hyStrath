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

#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "multi2Thermo.H" // NEW VINCENT
#include "addToRunTimeSelectionTable.H"
#include "fixed2TREnergyFvPatchScalarField.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fixed2TREnergyFvPatchScalarField::
fixed2TREnergyFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF)
{} // Only this constructor is used at run-time


Foam::fixed2TREnergyFvPatchScalarField::
fixed2TREnergyFvPatchScalarField
(
    const fixed2TREnergyFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper)
{}


Foam::fixed2TREnergyFvPatchScalarField::
fixed2TREnergyFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF, dict)
{}


Foam::fixed2TREnergyFvPatchScalarField::
fixed2TREnergyFvPatchScalarField
(
    const fixed2TREnergyFvPatchScalarField& tppsf
)
:
    fixedValueFvPatchScalarField(tppsf)
{}


Foam::fixed2TREnergyFvPatchScalarField::
fixed2TREnergyFvPatchScalarField
(
    const fixed2TREnergyFvPatchScalarField& tppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(tppsf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fixed2TREnergyFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    //Info << "fixed2TREnergy is used for patch called " << patch().name() << endl;

    const multi2Thermo& thermo = multi2Thermo::lookup2Thermo(*this);
    const label patchi = patch().index();

    const scalarField& pw = thermo.p().boundaryField()[patchi];

    fvPatchScalarField& Ttw =
        const_cast<fvPatchScalarField&>(thermo.T().boundaryField()[patchi]);
    Ttw.evaluate();

    operator==(thermo.het(pw, Ttw, patchi)); // Force an assignment, overriding fixedValue status

    fixedValueFvPatchScalarField::updateCoeffs();
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        fixed2TREnergyFvPatchScalarField
    );
}

// ************************************************************************* //
