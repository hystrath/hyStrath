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

#include "mixed2TREnergyFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "multi2Thermo.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::mixed2TREnergyFvPatchScalarField::
mixed2TREnergyFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(p, iF),
    thermo_(rho2ReactionThermo::lookup2ReactionThermo(*this))
{
    valueFraction() = 0.0;
    refValue() = 0.0;
    refGrad() = 0.0;
}


Foam::mixed2TREnergyFvPatchScalarField::
mixed2TREnergyFvPatchScalarField
(
    const mixed2TREnergyFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchScalarField(ptf, p, iF, mapper),
    thermo_(rho2ReactionThermo::lookup2ReactionThermo(*this))
{}


Foam::mixed2TREnergyFvPatchScalarField::
mixed2TREnergyFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchScalarField(p, iF, dict),
    thermo_(rho2ReactionThermo::lookup2ReactionThermo(*this))
{}


Foam::mixed2TREnergyFvPatchScalarField::
mixed2TREnergyFvPatchScalarField
(
    const mixed2TREnergyFvPatchScalarField& tppsf
)
:
    mixedFvPatchScalarField(tppsf),
    thermo_(rho2ReactionThermo::lookup2ReactionThermo(*this))
{}


Foam::mixed2TREnergyFvPatchScalarField::
mixed2TREnergyFvPatchScalarField
(
    const mixed2TREnergyFvPatchScalarField& tppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(tppsf, iF),
    thermo_(rho2ReactionThermo::lookup2ReactionThermo(*this))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::mixed2TREnergyFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const multi2Thermo& multiThermo = multi2Thermo::lookup2Thermo(*this);
    const label patchi = patch().index();

    const scalarField& pw = multiThermo.p().boundaryField()[patchi];
    mixedFvPatchScalarField& Tw = refCast<mixedFvPatchScalarField>
    (
        const_cast<fvPatchScalarField&>(multiThermo.T().boundaryField()[patchi])
    );
    Tw.evaluate();
    
    tmp<scalarField> thet(new scalarField(pw.size()));
    tmp<scalarField> thetRef(new scalarField(pw.size()));
    tmp<scalarField> thetfC(new scalarField(pw.size()));
    tmp<scalarField> tcvtTw(new scalarField(pw.size()));
    
    scalarField& het = thet.ref();
    scalarField& hetRef = thetRef.ref();
    scalarField& hetfC = thetfC.ref();
    scalarField& cvtTw = tcvtTw.ref();
    
    het = 0.0;
    hetRef = 0.0;
    hetfC = 0.0;
    cvtTw = 0.0;

    forAll(thermo_.composition().Y(), speciei)
    {
        fvPatchScalarField& spYw =
            const_cast<fvPatchScalarField&>
            (
                thermo_.composition().Y(speciei).boundaryField()[patchi]
            );
        spYw.evaluate();

        het += spYw*thermo_.composition().het(speciei, pw, Tw, patchi);
        hetRef += spYw
            *thermo_.composition().het(speciei, pw, Tw.refValue(), patchi);
        hetfC += spYw
            *thermo_.composition().het(speciei, pw, Tw, patch().faceCells());
        cvtTw += spYw*thermo_.composition().Cv_t(speciei, pw, Tw, patchi);
    }
    
    cvtTw *= Tw.refGrad();
    
    valueFraction() = Tw.valueFraction();
    refValue() = thetRef;
    refGrad() = tcvtTw + patch().deltaCoeffs()*(thet - thetfC);

    mixedFvPatchScalarField::updateCoeffs();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        mixed2TREnergyFvPatchScalarField
    );
}

// ************************************************************************* //
