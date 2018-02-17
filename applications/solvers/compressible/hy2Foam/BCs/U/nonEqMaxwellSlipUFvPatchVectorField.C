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

#include "nonEqMaxwellSlipUFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "mathematicalConstants.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "fvcGrad.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::nonEqMaxwellSlipUFvPatchVectorField::nonEqMaxwellSlipUFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    mixedFixedValueSlipFvPatchVectorField(p, iF),
    TName_("Tt"),
    rhoName_("rho"),
    muName_("mu"),
    tauMCName_("tauMC"),
    mfpName_("mfp"), // NEW VINCENT 28/02/2016
    accommodationCoeff_(1.0),
    vectorUwall_(vector::zero), // NEW VINCENT 21/10/2016
    //Uwall_(p.size(), vector::zero),
    Uwall_(p.size(), vectorUwall_), // NEW VINCENT 21/10/2016
    thermalCreep_(true),
    curvature_(true)
{}


Foam::nonEqMaxwellSlipUFvPatchVectorField::nonEqMaxwellSlipUFvPatchVectorField
(
    const nonEqMaxwellSlipUFvPatchVectorField& mspvf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFixedValueSlipFvPatchVectorField(mspvf, p, iF, mapper),
    TName_(mspvf.TName_),
    rhoName_(mspvf.rhoName_),
    muName_(mspvf.muName_),
    tauMCName_(mspvf.tauMCName_),
    mfpName_(mspvf.mfpName_), // NEW VINCENT 28/02/2016
    accommodationCoeff_(mspvf.accommodationCoeff_),
    vectorUwall_(mspvf.vectorUwall_), // NEW VINCENT 21/10/2016
    Uwall_(mspvf.Uwall_),
    thermalCreep_(mspvf.thermalCreep_),
    curvature_(mspvf.curvature_)
{}


Foam::nonEqMaxwellSlipUFvPatchVectorField::nonEqMaxwellSlipUFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFixedValueSlipFvPatchVectorField(p, iF),
    TName_(dict.lookupOrDefault<word>("Tt", "Tt")),
    rhoName_(dict.lookupOrDefault<word>("rho", "rho")),
    muName_(dict.lookupOrDefault<word>("mu", "mu")),
    tauMCName_(dict.lookupOrDefault<word>("tauMC", "tauMC")),
    mfpName_(dict.lookupOrDefault<word>("mfp", "mfp")), // NEW VINCENT 28/02/2016
    accommodationCoeff_(readScalar(dict.lookup("accommodationCoeff"))),
    vectorUwall_(dict.lookupOrDefault<vector>("Uwall", vector::zero)), // NEW VINCENT 21/10/2016
    //Uwall_("Uwall_", dict, p.size()),
    Uwall_(p.size(), vectorUwall_), // NEW VINCENT 21/10/2016
    thermalCreep_(dict.lookupOrDefault("thermalCreep", true)),
    curvature_(dict.lookupOrDefault("curvature", true))
{
    if
    (
        mag(accommodationCoeff_) < SMALL
     || mag(accommodationCoeff_) > 2.0
    )
    {
        FatalIOErrorIn
        (
            "maxwellSlipUFvPatchScalarField::maxwellSlipUFvPatchScalarField"
            "("
                "const fvPatch&, "
                "const DimensionedField<vector, volMesh>&, "
                "const dictionary&"
            ")",
            dict
        )   << "unphysical accommodationCoeff_ specified"
            << "(0 < accommodationCoeff_ <= 1)" << endl
            << exit(FatalIOError);
    }

    if (dict.found("value"))
    {
        fvPatchField<vector>::operator=
        (
            vectorField("value", dict, p.size())
        );

        if (dict.found("refValue") && dict.found("valueFraction"))
        {
            this->refValue() = vectorField("refValue", dict, p.size());
            this->valueFraction() =
                scalarField("valueFraction", dict, p.size());
        }
        else
        {
            this->refValue() = *this;
            this->valueFraction() = scalar(1.0);
        }
    }
    
}


Foam::nonEqMaxwellSlipUFvPatchVectorField::nonEqMaxwellSlipUFvPatchVectorField
(
    const nonEqMaxwellSlipUFvPatchVectorField& mspvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    mixedFixedValueSlipFvPatchVectorField(mspvf, iF),
    TName_(mspvf.TName_),
    rhoName_(mspvf.rhoName_),
    muName_(mspvf.muName_),
    tauMCName_(mspvf.tauMCName_),
    mfpName_(mspvf.mfpName_), // NEW VINCENT 28/02/2016
    accommodationCoeff_(mspvf.accommodationCoeff_),
    vectorUwall_(mspvf.vectorUwall_), // NEW VINCENT 21/10/2016
    Uwall_(mspvf.Uwall_),
    thermalCreep_(mspvf.thermalCreep_),
    curvature_(mspvf.curvature_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::nonEqMaxwellSlipUFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const fvPatchScalarField& pmu =
        patch().lookupPatchField<volScalarField, scalar>(muName_);
    const fvPatchScalarField& prho =
        patch().lookupPatchField<volScalarField, scalar>(rhoName_);
    const fvPatchScalarField& pmfp =
        patch().lookupPatchField<volScalarField, scalar>(mfpName_); // NEW VINCENT 28/02/2016    

    Field<scalar> pnu(pmu/prho);
    
    // DELETION VINCENT 28/02/2016 ********************************************
    /*Field<scalar> C1
    (
        sqrt(ppsi*constant::mathematical::piByTwo)
      * (2.0 - accommodationCoeff_)/accommodationCoeff_
    );*/
    // END DELETION VINCENT 28/02/2016 ****************************************
    Field<scalar> C1
    (
        pmfp/pnu
      * (2.0 - accommodationCoeff_)/accommodationCoeff_
    ); // NEW VINCENT 28/02/2016: divided by nu to keep the same units as before
    // Thus, what follows remains valid.
    
    valueFraction() = (1.0/(1.0 + patch().deltaCoeffs()*C1*pnu));

    refValue() = Uwall_;

    if (thermalCreep_)
    {
        const volScalarField& vsfT =
            this->db().objectRegistry::lookupObject<volScalarField>(TName_);
        label patchi = this->patch().index();
        const fvPatchScalarField& pT = vsfT.boundaryField()[patchi];
        Field<vector> gradpT(fvc::grad(vsfT)().boundaryField()[patchi]);
        vectorField n(patch().nf());

        refValue() -= 3.0*pnu/(4.0*pT)*transform(I - n*n, gradpT);
    }

    if (curvature_)
    {
        const fvPatchTensorField& ptauMC =
            patch().lookupPatchField<volTensorField, tensor>(tauMCName_);
        vectorField n(patch().nf());

        refValue() -= C1/prho*transform(I - n*n, (n & ptauMC));
    }

    mixedFixedValueSlipFvPatchVectorField::updateCoeffs();
}


void Foam::nonEqMaxwellSlipUFvPatchVectorField::write(Ostream& os) const
{
    fvPatchVectorField::write(os);
    writeEntryIfDifferent<word>(os, "Tt", "Tt", TName_);
    writeEntryIfDifferent<word>(os, "rho", "rho", rhoName_);
    writeEntryIfDifferent<word>(os, "mu", "mu", muName_);
    writeEntryIfDifferent<word>(os, "tauMC", "tauMC", tauMCName_);
    writeEntryIfDifferent<word>(os, "mfp", "mfp", mfpName_); // NEW VINCENT 28/02/2016

    os.writeKeyword("accommodationCoeff")
        << accommodationCoeff_ << token::END_STATEMENT << nl;
    os.writeKeyword("Uwall")
        << vectorUwall_ << token::END_STATEMENT << nl;  // NEW VINCENT 21/10/2016
    //Uwall_.writeEntry("Uwall", os); // DELETED VINCENT 21/10/2016
    os.writeKeyword("thermalCreep")
        << thermalCreep_ << token::END_STATEMENT << nl;
    os.writeKeyword("curvature") << curvature_ << token::END_STATEMENT << nl;
    /*os.writeKeyword("refValue")
        << refValue() << token::END_STATEMENT << nl;*/
    //refValue().writeEntry("refValue", os);
    //valueFraction().writeEntry("valueFraction", os); // DELETED VINCENT 21/10/2016

    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchVectorField,
        nonEqMaxwellSlipUFvPatchVectorField
    );
}

// ************************************************************************* //
