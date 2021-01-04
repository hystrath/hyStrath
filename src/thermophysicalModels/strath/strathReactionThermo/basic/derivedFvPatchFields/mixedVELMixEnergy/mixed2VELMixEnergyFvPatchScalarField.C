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

#include "mixed2VELMixEnergyFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "multi2Thermo.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::mixed2VELMixEnergyFvPatchScalarField::
mixed2VELMixEnergyFvPatchScalarField
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


Foam::mixed2VELMixEnergyFvPatchScalarField::
mixed2VELMixEnergyFvPatchScalarField
(
    const mixed2VELMixEnergyFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchScalarField(ptf, p, iF, mapper),
    thermo_(rho2ReactionThermo::lookup2ReactionThermo(*this))
{}


Foam::mixed2VELMixEnergyFvPatchScalarField::
mixed2VELMixEnergyFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchScalarField(p, iF, dict),
    thermo_(rho2ReactionThermo::lookup2ReactionThermo(*this))
{}


Foam::mixed2VELMixEnergyFvPatchScalarField::
mixed2VELMixEnergyFvPatchScalarField
(
    const mixed2VELMixEnergyFvPatchScalarField& tppsf
)
:
    mixedFvPatchScalarField(tppsf),
    thermo_(rho2ReactionThermo::lookup2ReactionThermo(*this))
{}


Foam::mixed2VELMixEnergyFvPatchScalarField::
mixed2VELMixEnergyFvPatchScalarField
(
    const mixed2VELMixEnergyFvPatchScalarField& tppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(tppsf, iF),
    thermo_(rho2ReactionThermo::lookup2ReactionThermo(*this))
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

    tmp<scalarField> thevel(new scalarField(pw.size()));
    tmp<scalarField> thevelRef(new scalarField(pw.size()));
    tmp<scalarField> thevelfC(new scalarField(pw.size()));
    tmp<scalarField> tcvvelTvw(new scalarField(pw.size()));
    
    scalarField& hevel = thevel.ref();
    scalarField& hevelRef = thevelRef.ref();
    scalarField& hevelfC = thevelfC.ref();
    scalarField& cvvelTvw = tcvvelTvw.ref();
    
    hevel = 0.0;
    hevelRef = 0.0;
    hevelfC = 0.0;
    cvvelTvw = 0.0;

    forAll(thermo_.composition().Y(), speciei)
    {
        fvPatchScalarField& spYw =
            const_cast<fvPatchScalarField&>
            (
                thermo_.composition().Y(speciei).boundaryField()[patchi]
            );
        spYw.evaluate();

        hevel += spYw*thermo_.composition().hevel(speciei, pw, Tvw, patchi);
        hevelRef += spYw
            *thermo_.composition().hevel(speciei, pw, Tvw.refValue(), patchi);
        hevelfC += spYw
            *thermo_.composition().hevel(speciei, pw, Tvw, patch().faceCells());
        cvvelTvw += spYw*thermo_.composition().Cv_vel(speciei, pw, Tvw, patchi);
    }
    
    cvvelTvw *= Tvw.refGrad();

    refValue() = thevelRef;
    refGrad() = tcvvelTvw + patch().deltaCoeffs()*(thevel - thevelfC);

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
