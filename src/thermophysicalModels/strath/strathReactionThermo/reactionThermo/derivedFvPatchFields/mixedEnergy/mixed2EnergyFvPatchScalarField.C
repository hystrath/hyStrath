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

#include "mixed2EnergyFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "multi2Thermo.H" // NEW VINCENT

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::mixed2EnergyFvPatchScalarField::
mixed2EnergyFvPatchScalarField
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


Foam::mixed2EnergyFvPatchScalarField::
mixed2EnergyFvPatchScalarField
(
    const mixed2EnergyFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchScalarField(ptf, p, iF, mapper),
    thermo_(rho2ReactionThermo::lookup2ReactionThermo(*this))
{}


Foam::mixed2EnergyFvPatchScalarField::
mixed2EnergyFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchScalarField(p, iF, dict),
    thermo_(rho2ReactionThermo::lookup2ReactionThermo(*this))
{}


Foam::mixed2EnergyFvPatchScalarField::
mixed2EnergyFvPatchScalarField
(
    const mixed2EnergyFvPatchScalarField& tppsf
)
:
    mixedFvPatchScalarField(tppsf),
    thermo_(rho2ReactionThermo::lookup2ReactionThermo(*this))
{}


Foam::mixed2EnergyFvPatchScalarField::
mixed2EnergyFvPatchScalarField
(
    const mixed2EnergyFvPatchScalarField& tppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(tppsf, iF),
    thermo_(rho2ReactionThermo::lookup2ReactionThermo(*this))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::mixed2EnergyFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const multi2Thermo& multiThermo = multi2Thermo::lookup2Thermo(*this);
    const label patchi = patch().index();

    const scalarField& pw = multiThermo.p().boundaryField()[patchi];

    mixedFvPatchScalarField& Ttw = refCast<mixedFvPatchScalarField>
    (
        const_cast<fvPatchScalarField&>(multiThermo.T().boundaryField()[patchi])
    );
    Ttw.evaluate();

    tmp<Field<scalar> > thevel(new Field<scalar>(pw.size()));
    Field<scalar>& hevel = thevel.ref();
    hevel = 0.0;

    tmp<Field<scalar> > thevelRef(new Field<scalar>(pw.size()));
    Field<scalar>& hevelRef = thevelRef.ref();
    hevelRef = 0.0;

    tmp<Field<scalar> > thevelfC(new Field<scalar>(pw.size()));
    Field<scalar>& hevelfC = thevelfC.ref();
    hevelfC = 0.0;

    tmp<Field<scalar> > tCvvelTvw(new Field<scalar>(pw.size()));
    Field<scalar>& cvvelTvw = tCvvelTvw.ref();
    cvvelTvw = 0.0;

    for(label speciei=0 ; speciei<thermo_.composition().Y().size() ; speciei++)
    {
        fvPatchScalarField& spYw =
            const_cast<fvPatchScalarField&>(thermo_.composition().Y(speciei).boundaryField()[patchi]);
        spYw.evaluate();

        mixedFvPatchScalarField& spTvw = refCast<mixedFvPatchScalarField>
            (const_cast<fvPatchScalarField&>(thermo_.composition().Tv(speciei).boundaryField()[patchi]));
        spTvw.evaluate();

        hevel += spYw*thermo_.composition().hevel(speciei, pw, spTvw, patchi);
        hevelRef += spYw*thermo_.composition().hevel(speciei, pw, spTvw.refValue(), patchi);
        hevelfC += spYw*thermo_.composition().hevel(speciei, pw, spTvw, patch().faceCells());
        cvvelTvw += spYw*thermo_.composition().Cv_vel(speciei, pw, spTvw, patchi)*spTvw.refGrad();
    }

    valueFraction() = Ttw.valueFraction();
    refValue() = multiThermo.het(pw, Ttw.refValue(), patchi) + thevelRef;
    refGrad() = multiThermo.Cv_t(pw, Ttw, patchi)*Ttw.refGrad() + tCvvelTvw
        + patch().deltaCoeffs()*
          (
              multiThermo.het(pw, Ttw, patchi) + thevel
            - multiThermo.het(pw, Ttw, patch().faceCells()) - thevelfC
          );

    mixedFvPatchScalarField::updateCoeffs();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        mixed2EnergyFvPatchScalarField
    );
}

// ************************************************************************* //
