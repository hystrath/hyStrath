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

#include "lumpedMassWallTemperatureFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "mappedPatchBase.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::lumpedMassWallTemperatureFvPatchScalarField::
lumpedMassWallTemperatureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(p, iF),
    temperatureCoupledBase(patch(), "undefined", "undefined", "undefined-K"),
    Cp_(0.0),
    mass_(0.0),
    curTimeIndex_(-1)
{
    refValue() = 0.0;
    refGrad() = 0.0;
    valueFraction() = 1.0;
}


Foam::lumpedMassWallTemperatureFvPatchScalarField::
lumpedMassWallTemperatureFvPatchScalarField
(
    const lumpedMassWallTemperatureFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchScalarField(ptf, p, iF, mapper),
    temperatureCoupledBase(patch(), ptf),
    Cp_(ptf.Cp_),
    mass_(ptf.mass_),
    curTimeIndex_(-1)
{}


Foam::lumpedMassWallTemperatureFvPatchScalarField::
lumpedMassWallTemperatureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchScalarField(p, iF),
    temperatureCoupledBase(patch(), dict),
    Cp_(readScalar(dict.lookup("Cp"))),
    mass_(readScalar(dict.lookup("mass"))),
    curTimeIndex_(-1)
{
    refGrad() = 0.0;
    valueFraction() = 1.0;
    refValue() = scalarField("value", dict, p.size());

    fvPatchScalarField::operator=(scalarField("value", dict, p.size()));
}


Foam::lumpedMassWallTemperatureFvPatchScalarField::
lumpedMassWallTemperatureFvPatchScalarField
(
    const lumpedMassWallTemperatureFvPatchScalarField& tppsf
)
:
    mixedFvPatchScalarField(tppsf),
    temperatureCoupledBase(tppsf),
    Cp_(tppsf.Cp_),
    mass_(tppsf.mass_),
    curTimeIndex_(-1)
{}


Foam::lumpedMassWallTemperatureFvPatchScalarField::
lumpedMassWallTemperatureFvPatchScalarField
(
    const lumpedMassWallTemperatureFvPatchScalarField& tppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(tppsf, iF),
    temperatureCoupledBase(patch(), tppsf),
    Cp_(tppsf.Cp_),
    mass_(tppsf.mass_),
    curTimeIndex_(-1)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::lumpedMassWallTemperatureFvPatchScalarField::updateCoeffs()
{
    if (updated() || (curTimeIndex_ == this->db().time().timeIndex()))
    {
        return;
    }

    // Calculate heat flux in or out the wall
    scalarField& Tp(*this);
    const scalarField& magSf = patch().magSf();

    const scalar deltaT(db().time().deltaTValue());

    tmp<scalarField> tkappa(kappa(Tp));

    const scalarField q(tkappa.ref()*snGrad());

    // Total heat in or out of the wall
    const scalar Q = gSum(q*magSf);

    Tp += -(Q/mass_/Cp_)*deltaT;

    refGrad() = 0.0;
    refValue() = Tp;
    valueFraction() = 1.0;

    mixedFvPatchScalarField::updateCoeffs();

    if (debug)
    {
        scalar qin(0);
        scalar qout(0);

        forAll(q, facei)
        {
            if (q[facei] > 0.0) //out the wall
            {
                qout += q[facei]*magSf[facei];
            }
            else if (q[facei] < 0.0) //into the wall
            {
                qin += q[facei]*magSf[facei];
            }
        }

        Info<< patch().boundaryMesh().mesh().name() << ':'
            << patch().name() << ':'
            << this->internalField().name() << " :"
            << " heat transfer rate:" << Q
            << " wall temperature "
            << " min:" << gMin(*this)
            << " max:" << gMax(*this)
            << " avg:" << gAverage(*this)
            << " Qin [W]:" << qin
            << " Qout [W]:" << qout
            << endl;
    }

    curTimeIndex_ = this->db().time().timeIndex();
}


void Foam::lumpedMassWallTemperatureFvPatchScalarField::write
(
    Ostream& os
) const
{
    mixedFvPatchScalarField::write(os);
    temperatureCoupledBase::write(os);

    os.writeKeyword("Cp")<< Cp_ << token::END_STATEMENT << nl;
    os.writeKeyword("mass")<< mass_
        << token::END_STATEMENT << nl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        lumpedMassWallTemperatureFvPatchScalarField
    );
}

// ************************************************************************* //
