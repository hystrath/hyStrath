/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012-2014 OpenFOAM Foundation
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

#include "omegaLowReWallFunctionFvPatchScalarField.H"
#include "strath/compressible/turbulenceModel/turbulenceModel/turbulenceModel2.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"
#include "wallFvPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void omegaLowReWallFunctionFvPatchScalarField::checkType()
{
    if (!isA<wallFvPatch>(patch()))
    {
        FatalErrorIn("omegaLowReWallFunctionFvPatchScalarField::checkType()")
            << "Invalid wall function specification" << nl
            << "    Patch type for patch " << patch().name()
            << " must be wall" << nl
            << "    Current patch type is " << patch().type() << nl << endl
            << abort(FatalError);
    }
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

scalar omegaLowReWallFunctionFvPatchScalarField::yPlusLam
(
    const scalar kappa,
    const scalar E
)
{
    scalar ypl = 11.0;

    for (int i=0; i<10; i++)
    {
        ypl = log(max(E*ypl, 1))/kappa;
    }

    return ypl;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

omegaLowReWallFunctionFvPatchScalarField::omegaLowReWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchField<scalar>(p, iF),
    beta1_(0.09),
    kappa_(0.41),
    E_(9.8),
    yPlusLam_(yPlusLam(kappa_, E_))
{
    checkType();
}


omegaLowReWallFunctionFvPatchScalarField::omegaLowReWallFunctionFvPatchScalarField
(
    const omegaLowReWallFunctionFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<scalar>(ptf, p, iF, mapper),
    beta1_(ptf.beta1_),
    kappa_(ptf.kappa_),
    E_(ptf.E_),
    yPlusLam_(ptf.yPlusLam_)
{
    checkType();
}


omegaLowReWallFunctionFvPatchScalarField::omegaLowReWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<scalar>(p, iF, dict),
    beta1_(dict.lookupOrDefault<scalar>("beta1", 0.075)),
    kappa_(dict.lookupOrDefault<scalar>("kappa", 0.41)),
    E_(dict.lookupOrDefault<scalar>("E", 9.8)),
    yPlusLam_(yPlusLam(kappa_, E_))
{
    checkType();
}


omegaLowReWallFunctionFvPatchScalarField::omegaLowReWallFunctionFvPatchScalarField
(
    const omegaLowReWallFunctionFvPatchScalarField& owfpsf
)
:
    fixedValueFvPatchField<scalar>(owfpsf),
    beta1_(owfpsf.beta1_),
    kappa_(owfpsf.kappa_),
    E_(owfpsf.E_),
    yPlusLam_(owfpsf.yPlusLam_)
{
    checkType();
}


omegaLowReWallFunctionFvPatchScalarField::omegaLowReWallFunctionFvPatchScalarField
(
    const omegaLowReWallFunctionFvPatchScalarField& owfpsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchField<scalar>(owfpsf, iF),
    beta1_(owfpsf.beta1_),
    kappa_(owfpsf.kappa_),
    E_(owfpsf.E_),
    yPlusLam_(owfpsf.yPlusLam_)
{
    checkType();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void omegaLowReWallFunctionFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const label patchI = patch().index();

    const turbulenceModel2& turbulence =
        db().lookupObject<turbulenceModel2>("turbulenceModel2");
    const scalarField& y = turbulence.y()[patchI];

    scalarField& ow = *this;

    const tmp<volScalarField> tmu = turbulence.mu();
    const scalarField& muw = tmu().boundaryField()[patchI];

    const scalarField& rhow = turbulence.rho().boundaryField()[patchI];

    // Set omega wall values
    forAll(ow, faceI)
    {
        scalar nuw = muw[faceI]/rhow[faceI];

        ow[faceI] = 60.0*nuw/(beta1_*pow(y[faceI], 2.0));
    }

    fixedValueFvPatchField<scalar>::updateCoeffs();

    // TODO: perform averaging for cells sharing more than one boundary face
}


void omegaLowReWallFunctionFvPatchScalarField::evaluate
(
    const Pstream::commsTypes commsType
)
{
    fixedValueFvPatchField<scalar>::evaluate(commsType);
}


void omegaLowReWallFunctionFvPatchScalarField::write(Ostream& os) const
{
    fvPatchField<scalar>::write(os);
    os.writeKeyword("beta1") << beta1_ << token::END_STATEMENT << nl;
    os.writeKeyword("kappa") << kappa_ << token::END_STATEMENT << nl;
    os.writeKeyword("E") << E_ << token::END_STATEMENT << nl;
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    omegaLowReWallFunctionFvPatchScalarField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace compressible
} // End namespace Foam

// ************************************************************************* //
