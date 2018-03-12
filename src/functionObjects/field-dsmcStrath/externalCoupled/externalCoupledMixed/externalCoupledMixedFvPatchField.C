/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013-2015 OpenFOAM Foundation
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

#include "externalCoupledMixedFvPatchField.H"
#include "fvPatchFieldMapper.H"
#include "ISstream.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::externalCoupledMixedFvPatchField<Type>::
externalCoupledMixedFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    mixedFvPatchField<Type>(p, iF)
{
    this->refValue() = Type(Zero);
    this->refGrad() = Type(Zero);
    this->valueFraction() = 0.0;
}


template<class Type>
Foam::externalCoupledMixedFvPatchField<Type>::
externalCoupledMixedFvPatchField
(
    const externalCoupledMixedFvPatchField& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchField<Type>(ptf, p, iF, mapper)
{}


template<class Type>
Foam::externalCoupledMixedFvPatchField<Type>::
externalCoupledMixedFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchField<Type>(p, iF, dict)
{}


template<class Type>
Foam::externalCoupledMixedFvPatchField<Type>::
externalCoupledMixedFvPatchField
(
    const externalCoupledMixedFvPatchField& ecmpf
)
:
    mixedFvPatchField<Type>(ecmpf)
{}


template<class Type>
Foam::externalCoupledMixedFvPatchField<Type>::
externalCoupledMixedFvPatchField
(
    const externalCoupledMixedFvPatchField& ecmpf,
    const DimensionedField<Type, volMesh>& iF
)
:
    mixedFvPatchField<Type>(ecmpf, iF)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
Foam::externalCoupledMixedFvPatchField<Type>::
~externalCoupledMixedFvPatchField()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::externalCoupledMixedFvPatchField<Type>::writeHeader
(
    Ostream& os
) const
{
    os  << "# Values: value snGrad refValue refGrad valueFraction" << endl;
}


template<class Type>
void Foam::externalCoupledMixedFvPatchField<Type>::writeData
(
    Ostream& os
) const
{
    const Field<Type> snGrad(this->snGrad());
    const Field<Type>& refValue(this->refValue());
    const Field<Type>& refGrad(this->refGrad());
    const scalarField& valueFraction(this->valueFraction());

    forAll(refValue, facei)
    {
        os  << this->operator[](facei) << token::SPACE
            << snGrad[facei] << token::SPACE
            << refValue[facei] << token::SPACE
            << refGrad[facei] << token::SPACE
            << valueFraction[facei] << nl;
    }
}


template<class Type>
void Foam::externalCoupledMixedFvPatchField<Type>::readData(Istream& is)
{
    // Assume generic input stream so we can do line-based format and skip
    // unused columns
    ISstream& iss = dynamic_cast<ISstream&>(is);

    string line;

    forAll(*this, facei)
    {
        iss.getLine(line);
        IStringStream lineStr(line);

        // For symmetry with writing ignore value, snGrad columns

        Type value, snGrad;

        lineStr
            >> value
            >> snGrad
            >> this->refValue()[facei]
            >> this->refGrad()[facei]
            >> this->valueFraction()[facei];
    }
}


// ************************************************************************* //
