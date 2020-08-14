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

#include "gradient2EnergyFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "multi2Thermo.H" // NEW VINCENT

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::gradient2EnergyFvPatchScalarField::
gradient2EnergyFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(p, iF),
    thermo_(rho2ReactionThermo::lookup2ReactionThermo(*this))
{} // Only this constructor is used at run-time


Foam::gradient2EnergyFvPatchScalarField::
gradient2EnergyFvPatchScalarField
(
    const gradient2EnergyFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedGradientFvPatchScalarField(ptf, p, iF, mapper),
    thermo_(rho2ReactionThermo::lookup2ReactionThermo(*this))
{}


Foam::gradient2EnergyFvPatchScalarField::
gradient2EnergyFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedGradientFvPatchScalarField(p, iF, dict),
    thermo_(rho2ReactionThermo::lookup2ReactionThermo(*this))
{}


Foam::gradient2EnergyFvPatchScalarField::
gradient2EnergyFvPatchScalarField
(
    const gradient2EnergyFvPatchScalarField& tppsf
)
:
    fixedGradientFvPatchScalarField(tppsf),
    thermo_(rho2ReactionThermo::lookup2ReactionThermo(*this))
{}


Foam::gradient2EnergyFvPatchScalarField::
gradient2EnergyFvPatchScalarField
(
    const gradient2EnergyFvPatchScalarField& tppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(tppsf, iF),
    thermo_(rho2ReactionThermo::lookup2ReactionThermo(*this))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::gradient2EnergyFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    //Info << "gradient2Energy is used for patch called " << patch().name() << endl;

    const multi2Thermo& multiThermo = multi2Thermo::lookup2Thermo(*this);
    const label patchi = patch().index();

    const scalarField& pw = multiThermo.p().boundaryField()[patchi];

    fvPatchScalarField& Ttw =
        const_cast<fvPatchScalarField&>(multiThermo.T().boundaryField()[patchi]);
    Ttw.evaluate();

    tmp<Field<scalar> > thevel(new Field<scalar>(pw.size()));
    Field<scalar>& hevel = thevel.ref();
    hevel = 0.0;

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

        fvPatchScalarField& spTvw =
            const_cast<fvPatchScalarField&>(thermo_.composition().Tv(speciei).boundaryField()[patchi]);
        spTvw.evaluate();

        hevel += spYw*thermo_.composition().hevel(speciei, pw, spTvw, patchi);
        hevelfC += spYw*thermo_.composition().hevel(speciei, pw, spTvw, patch().faceCells());
        cvvelTvw += spYw*thermo_.composition().Cv_vel(speciei, pw, spTvw, patchi)*spTvw.snGrad();
    }

    gradient() = multiThermo.Cv_t(pw, Ttw, patchi)*Ttw.snGrad() + tCvvelTvw
      + patch().deltaCoeffs()*
        (
            multiThermo.het(pw, Ttw, patchi) + thevel
          - multiThermo.het(pw, Ttw, patch().faceCells()) - thevelfC
        );

    fixedGradientFvPatchScalarField::updateCoeffs();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        gradient2EnergyFvPatchScalarField
    );
}

// ************************************************************************* //
