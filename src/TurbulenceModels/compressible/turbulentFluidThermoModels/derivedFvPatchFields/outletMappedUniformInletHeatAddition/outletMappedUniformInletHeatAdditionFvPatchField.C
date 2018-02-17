/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016 OpenCFD Ltd.
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

#include "outletMappedUniformInletHeatAdditionFvPatchField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "basic2Thermo.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::outletMappedUniformInletHeatAdditionFvPatchField::
outletMappedUniformInletHeatAdditionFvPatchField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF),
    outletPatchName_(),
    phiName_("phi"),
    Q_(0),
    TMin_(0),
    TMax_(5000)
{}


Foam::outletMappedUniformInletHeatAdditionFvPatchField::
outletMappedUniformInletHeatAdditionFvPatchField
(
    const outletMappedUniformInletHeatAdditionFvPatchField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper),
    outletPatchName_(ptf.outletPatchName_),
    phiName_(ptf.phiName_),
    Q_(ptf.Q_),
    TMin_(ptf.TMin_),
    TMax_(ptf.TMax_)
{}


Foam::outletMappedUniformInletHeatAdditionFvPatchField::
outletMappedUniformInletHeatAdditionFvPatchField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF, dict),
    outletPatchName_(dict.lookup("outletPatch")),
    phiName_(dict.lookupOrDefault<word>("phi", "phi")),
    Q_(readScalar(dict.lookup("Q"))),
    TMin_(dict.lookupOrDefault<scalar>("TMin", 0)),
    TMax_(dict.lookupOrDefault<scalar>("TMax", 5000))
{}



Foam::outletMappedUniformInletHeatAdditionFvPatchField::
outletMappedUniformInletHeatAdditionFvPatchField
(
    const outletMappedUniformInletHeatAdditionFvPatchField& ptf
)
:
    fixedValueFvPatchScalarField(ptf),
    outletPatchName_(ptf.outletPatchName_),
    phiName_(ptf.phiName_),
    Q_(ptf.Q_),
    TMin_(ptf.TMin_),
    TMax_(ptf.TMax_)
{}



Foam::outletMappedUniformInletHeatAdditionFvPatchField::
outletMappedUniformInletHeatAdditionFvPatchField
(
    const outletMappedUniformInletHeatAdditionFvPatchField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(ptf, iF),
    outletPatchName_(ptf.outletPatchName_),
    phiName_(ptf.phiName_),
    Q_(ptf.Q_),
    TMin_(ptf.TMin_),
    TMax_(ptf.TMax_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void Foam::outletMappedUniformInletHeatAdditionFvPatchField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    const volScalarField& vsf =
    (
        dynamic_cast<const volScalarField&>(this->internalField())
    );

    const fvPatch& fvp = this->patch();

    label outletPatchID =
        fvp.patch().boundaryMesh().findPatchID(outletPatchName_);

    if (outletPatchID < 0)
    {
        FatalErrorInFunction
            << "Unable to find outlet patch " << outletPatchName_
            << abort(FatalError);
    }

    const fvPatch& outletPatch = fvp.boundaryMesh()[outletPatchID];

    const fvPatchField<scalar>& outletPatchField =
        vsf.boundaryField()[outletPatchID];

    const surfaceScalarField& phi =
        db().lookupObject<surfaceScalarField>(phiName_);

    const scalarField& outletPatchPhi = phi.boundaryField()[outletPatchID];
    scalar sumOutletPatchPhi = gSum(outletPatchPhi);

    if (sumOutletPatchPhi > SMALL)
    {
        const basic2Thermo& thermo =
            db().lookupObject<basic2Thermo>(basic2Thermo::dictName);

        const scalarField& pp = thermo.p().boundaryField()[outletPatchID];
        const scalarField& pT = thermo.T().boundaryField()[outletPatchID];

        scalar averageOutletField =
            gSum(outletPatchPhi*outletPatchField)/sumOutletPatchPhi;

        const scalarField Cpf(thermo.Cp(pp, pT, outletPatchID));

        scalar totalPhiCp = gSum(outletPatchPhi)*gAverage(Cpf);

        operator==(min(max(averageOutletField + Q_/totalPhiCp, TMin_), TMax_));
    }
    else
    {
        scalar averageOutletField =
            gSum(outletPatch.magSf()*outletPatchField)
           /gSum(outletPatch.magSf());

        operator==(averageOutletField);
    }

    fixedValueFvPatchScalarField::updateCoeffs();
}


void Foam::outletMappedUniformInletHeatAdditionFvPatchField::write
(
    Ostream& os
) const
{
    fvPatchScalarField::write(os);
    os.writeKeyword("outletPatch")
        << outletPatchName_ << token::END_STATEMENT << nl;
    writeEntryIfDifferent<word>(os, "phi", "phi", phiName_);
    os.writeKeyword("Q") << Q_ << token::END_STATEMENT << nl;
    os.writeKeyword("TMin") << TMin_ << token::END_STATEMENT << nl;
    os.writeKeyword("TMax") << TMax_ << token::END_STATEMENT << nl;

    this->writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        outletMappedUniformInletHeatAdditionFvPatchField
    );
}


// ************************************************************************* //
