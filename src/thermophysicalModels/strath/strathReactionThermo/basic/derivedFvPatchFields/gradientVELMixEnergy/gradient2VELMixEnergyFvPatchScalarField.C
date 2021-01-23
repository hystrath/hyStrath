/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016-2021 hyStrath
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

#include "gradient2VELMixEnergyFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "multi2Thermo.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::gradient2VELMixEnergyFvPatchScalarField::
gradient2VELMixEnergyFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(p, iF),
    thermo_(rho2ReactionThermo::lookup2ReactionThermo(*this))
{} // Only this constructor is used


Foam::gradient2VELMixEnergyFvPatchScalarField::
gradient2VELMixEnergyFvPatchScalarField
(
    const gradient2VELMixEnergyFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedGradientFvPatchScalarField(ptf, p, iF, mapper),
    thermo_(rho2ReactionThermo::lookup2ReactionThermo(*this))
{}


Foam::gradient2VELMixEnergyFvPatchScalarField::
gradient2VELMixEnergyFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedGradientFvPatchScalarField(p, iF, dict),
    thermo_(rho2ReactionThermo::lookup2ReactionThermo(*this))
{}


Foam::gradient2VELMixEnergyFvPatchScalarField::
gradient2VELMixEnergyFvPatchScalarField
(
    const gradient2VELMixEnergyFvPatchScalarField& tppsf
)
:
    fixedGradientFvPatchScalarField(tppsf),
    thermo_(rho2ReactionThermo::lookup2ReactionThermo(*this))
{}


Foam::gradient2VELMixEnergyFvPatchScalarField::
gradient2VELMixEnergyFvPatchScalarField
(
    const gradient2VELMixEnergyFvPatchScalarField& tppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(tppsf, iF),
    thermo_(rho2ReactionThermo::lookup2ReactionThermo(*this))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::gradient2VELMixEnergyFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    //Info<< "gradient2VELMixEnergy is used for patch called "
    //    << patch().name() << endl;

    const multi2Thermo& multiThermo = multi2Thermo::lookup2Thermo(*this);
    const label patchi = patch().index();

    const scalarField& pw = multiThermo.p().boundaryField()[patchi];

    fvPatchScalarField& Tvw =
        const_cast<fvPatchScalarField&>
        (
            multiThermo.Tv().boundaryField()[patchi]
        );
    Tvw.evaluate();

    tmp<scalarField> thevel(new scalarField(pw.size()));
    tmp<scalarField> thevelfC(new scalarField(pw.size()));
    tmp<scalarField> tCvvelTvw(new scalarField(pw.size()));
    
    scalarField& hevel = thevel.ref();
    scalarField& hevelfC = thevelfC.ref();
    scalarField& cvvelTvw = tCvvelTvw.ref();
    
    hevel = 0.0;
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

        hevel += spYw
            *thermo_.composition().hevel(speciei, pw, Tvw, patchi);
        hevelfC += spYw
            *thermo_.composition().hevel(speciei, pw, Tvw, patch().faceCells());
        cvvelTvw += spYw
            *thermo_.composition().Cv_vel(speciei, pw, Tvw, patchi);
    }
    
    cvvelTvw *= Tvw.snGrad();

    gradient() = tCvvelTvw + patch().deltaCoeffs()*(thevel - thevelfC);

    fixedGradientFvPatchScalarField::updateCoeffs();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        gradient2VELMixEnergyFvPatchScalarField
    );
}

// ************************************************************************* //
