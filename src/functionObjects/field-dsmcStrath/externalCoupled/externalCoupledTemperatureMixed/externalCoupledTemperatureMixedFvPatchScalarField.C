/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013-2016 OpenFOAM Foundation
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

#include "externalCoupledTemperatureMixedFvPatchScalarField.H"
#include "turbulentFluidThermoModel.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::externalCoupledTemperatureMixedFvPatchScalarField::writeHeader
(
    Ostream& os
) const
{
    os  << "# Values: magSf T qDot htc" << endl;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::externalCoupledTemperatureMixedFvPatchScalarField::
externalCoupledTemperatureMixedFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    externalCoupledMixedFvPatchField<scalar>(p, iF)
{}


Foam::externalCoupledTemperatureMixedFvPatchScalarField::
externalCoupledTemperatureMixedFvPatchScalarField
(
    const externalCoupledTemperatureMixedFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    externalCoupledMixedFvPatchField<scalar>(ptf, p, iF, mapper)
{}


Foam::externalCoupledTemperatureMixedFvPatchScalarField::
externalCoupledTemperatureMixedFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    //externalCoupledMixedFvPatchField<scalar>(p, iF, dict)
    externalCoupledMixedFvPatchField<scalar>(p, iF)
{
    if (dict.found("refValue"))
    {
        // Initialise same way as mixed
        this->refValue() = scalarField("refValue", dict, p.size());
        this->refGrad() = scalarField("refGradient", dict, p.size());
        this->valueFraction() = scalarField("valueFraction", dict, p.size());

        evaluate();
    }
    else
    {
        // For convenience: initialise as fixedValue with either read value
        // or extrapolated value
        if (dict.found("value"))
        {
            fvPatchField<scalar>::operator=
            (
                scalarField("value", dict, p.size())
            );
        }
        else
        {
            fvPatchField<scalar>::operator=(this->patchInternalField());
        }

        // Initialise as a fixed value
        this->refValue() = *this;
        this->refGrad() = Zero;
        this->valueFraction() = 1.0;
    }
}


Foam::externalCoupledTemperatureMixedFvPatchScalarField::
externalCoupledTemperatureMixedFvPatchScalarField
(
    const externalCoupledTemperatureMixedFvPatchScalarField& ecmpf
)
:
    externalCoupledMixedFvPatchField<scalar>(ecmpf)
{}


Foam::externalCoupledTemperatureMixedFvPatchScalarField::
externalCoupledTemperatureMixedFvPatchScalarField
(
    const externalCoupledTemperatureMixedFvPatchScalarField& ecmpf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    externalCoupledMixedFvPatchField<scalar>(ecmpf, iF)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::externalCoupledTemperatureMixedFvPatchScalarField::
~externalCoupledTemperatureMixedFvPatchScalarField()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::externalCoupledTemperatureMixedFvPatchScalarField::writeData
(
    Ostream& os
) const
{
    const label patchi = patch().index();

    // Heat flux [W/m2]
    scalarField qDot(this->patch().size(), 0.0);

    typedef compressible::turbulenceModel cmpTurbModelType;

    static word turbName
    (
        IOobject::groupName
        (
            turbulenceModel::propertiesName,
            internalField().group()
        )
    );

    static word thermoName("thermophysicalProperties");

    if (db().foundObject<cmpTurbModelType>(turbName))
    {
        const cmpTurbModelType& turbModel =
            db().lookupObject<cmpTurbModelType>(turbName);

        const basicThermo& thermo = turbModel.transport();

        const fvPatchScalarField& hep = thermo.he().boundaryField()[patchi];

        qDot = turbModel.alphaEff(patchi)*hep.snGrad();
    }
    else if (db().foundObject<basicThermo>(thermoName))
    {
        const basicThermo& thermo = db().lookupObject<basicThermo>(thermoName);

        const fvPatchScalarField& hep = thermo.he().boundaryField()[patchi];

        qDot = thermo.alpha().boundaryField()[patchi]*hep.snGrad();
    }
    else
    {
        FatalErrorInFunction
            << "Condition requires either compressible turbulence and/or "
            << "thermo model to be available" << exit(FatalError);
    }

    // Patch temperature [K]
    const scalarField& Tp(*this);

    // Near wall cell temperature [K]
    const scalarField Tc(patchInternalField());

    // Heat transfer coefficient [W/m2/K]
    const scalarField htc(qDot/(Tp - Tc + ROOTVSMALL));

    const Field<scalar>& magSf(this->patch().magSf());

    forAll(patch(), facei)
    {
        os  << magSf[facei] << token::SPACE
            << Tp[facei] << token::SPACE
            << qDot[facei] << token::SPACE
            << htc[facei] << token::SPACE
            << nl;
    }
}


void Foam::externalCoupledTemperatureMixedFvPatchScalarField::readData
(
    Istream& is
)
{
    // Assume generic input stream so we can do line-based format and skip
    // unused columns
    ISstream& iss = dynamic_cast<ISstream&>(is);

    string line;

    forAll(*this, facei)
    {
        iss.getLine(line);
        IStringStream lineStr(line);

        lineStr
            >> this->refValue()[facei]
            >> this->refGrad()[facei]
            >> this->valueFraction()[facei];
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        externalCoupledTemperatureMixedFvPatchScalarField
    );
}


// ************************************************************************* //
