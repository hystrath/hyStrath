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

#include "fvInsulatorPatch.H"
#include "constants.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fvInsulatorPatch::fvInsulatorPatch
(
    const fvPatch& patch,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(patch, iF),
    phiEName_("phiE"),
    wallQName_("wallQ"),
    dz_(0.0),
    epsilonM_(0.0),
    phiW_(patch.size(), 0.0)
{}


Foam::fvInsulatorPatch::fvInsulatorPatch
(
    const fvInsulatorPatch& ptf,
    const fvPatch& patch,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, patch, iF, mapper),
    phiEName_(ptf.phiEName_),
    wallQName_(ptf.wallQName_),
    dz_(ptf.dz_),
    epsilonM_(ptf.epsilonM_),
    phiW_(ptf.phiW_, mapper)
{}


Foam::fvInsulatorPatch::fvInsulatorPatch
(
    const fvPatch& patch,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(patch, iF),
    phiEName_(dict.lookupOrDefault<word>("phiE", "phiE")),
    wallQName_(dict.lookupOrDefault<word>("wallQ", "wallQ")),
    dz_(readScalar(dict.lookup("dz"))),
    epsilonM_(readScalar(dict.lookup("epsilonM"))),
    phiW_("phiW", dict, patch.size())
{
    if (dict.found("value"))
    {
        fvPatchField<scalar>::operator=
        (
            scalarField("value", dict, patch.size())
        );
    }
    else
    {
        fvPatchField<scalar>::operator=(phiW_);
    }
}


Foam::fvInsulatorPatch::fvInsulatorPatch
(
    const fvInsulatorPatch& tppsf
)
:
    fixedValueFvPatchScalarField(tppsf),
    phiEName_(tppsf.phiEName_),
    wallQName_(tppsf.wallQName_),
    dz_(tppsf.dz_),
    epsilonM_(tppsf.epsilonM_),
    phiW_(tppsf.phiW_)
{}


Foam::fvInsulatorPatch::fvInsulatorPatch
(
    const fvInsulatorPatch& tppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(tppsf, iF),
    phiEName_(tppsf.phiEName_),
    wallQName_(tppsf.wallQName_),
    dz_(tppsf.dz_),
    epsilonM_(tppsf.epsilonM_),
    phiW_(tppsf.phiW_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fvInsulatorPatch::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedValueFvPatchScalarField::autoMap(m);
    phiW_.autoMap(m);
}


void Foam::fvInsulatorPatch::rmap
(
    const fvPatchScalarField& ptf,
    const labelList& addr
)
{
    fixedValueFvPatchScalarField::rmap(ptf, addr);

    const fvInsulatorPatch& tiptf =
        refCast<const fvInsulatorPatch>(ptf);

    phiW_.rmap(tiptf.phiW_, addr);
}


void Foam::fvInsulatorPatch::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const fvPatchScalarField& wallQp =
        patch().lookupPatchField<volScalarField, scalar>(wallQName_);

    const scalarField& faceArea = patch().magSf();


    scalar eps0 = constant::electromagnetic::epsilon0.value(); // permitivity of space

    //- determine potential based on parallel plate capacitor model
    forAll(wallQp,pI)
    {
        scalar C = eps0*epsilonM_*faceArea[pI]/dz_;

        phiW_[pI] = wallQp[pI]/C;
    }

    Info << "           max(phiW)        = " << gMax(phiW_) << endl;
    Info << "           min(phiW)        = " << gMin(phiW_) << endl;
    //- update boundary based on new charge density distribution

    operator==(phiW_);

    fixedValueFvPatchScalarField::updateCoeffs();
}


void Foam::fvInsulatorPatch::write(Ostream& os) const
{
    fvPatchScalarField::write(os);
    writeEntryIfDifferent<word>(os, "phiE", "phiE", phiEName_);
    writeEntryIfDifferent<word>(os, "wallQ", "wallQ", wallQName_);
    os.writeKeyword("dz") << dz_ << token::END_STATEMENT << nl;
    os.writeKeyword("epsilonM") << epsilonM_ << token::END_STATEMENT << nl;
    phiW_.writeEntry("phiW", os);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        fvInsulatorPatch
    );
}

// ************************************************************************* //
