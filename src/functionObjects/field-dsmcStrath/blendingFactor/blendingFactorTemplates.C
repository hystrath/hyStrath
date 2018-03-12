/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013-2016 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2016 OpenCFD Ltd.
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

#include "gaussConvectionScheme.H"
#include "boundedConvectionScheme.H"
#include "blendedSchemeBase.H"
#include "fvcCellReduce.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
void Foam::functionObjects::blendingFactor::calcBlendingFactor
(
    const GeometricField<Type, fvPatchField, volMesh>& field,
    const typename fv::convectionScheme<Type>& cs
)
{
    if (!isA<fv::gaussConvectionScheme<Type>>(cs))
    {
        WarningInFunction
            << "Scheme for field " << field.name() << " is not a "
            << fv::gaussConvectionScheme<Type>::typeName
            << " scheme. Not calculating " << resultName_ << endl;

        return;
    }

    const fv::gaussConvectionScheme<Type>& gcs =
        refCast<const fv::gaussConvectionScheme<Type>>(cs);


    const surfaceInterpolationScheme<Type>& interpScheme = gcs.interpScheme();

    if (!isA<blendedSchemeBase<Type>>(interpScheme))
    {
        WarningInFunction
            << interpScheme.type() << " is not a blended scheme"
            << ". Not calculating " << resultName_ << endl;

        return;
    }

    // Retrieve the face-based blending factor
    const blendedSchemeBase<Type>& blendedScheme =
        refCast<const blendedSchemeBase<Type>>(interpScheme);
    const surfaceScalarField factorf(blendedScheme.blendingFactor(field));

    // Convert into vol field whose values represent the local face minima
    // Note:
    // - factor applied to 1st scheme, and (1-factor) to 2nd scheme
    // - not using the store(...) mechanism due to need to correct BCs
    volScalarField& indicator =
        const_cast<volScalarField&>
        (
            lookupObject<volScalarField>(resultName_)
        );
    indicator = 1 - fvc::cellReduce(factorf, minEqOp<scalar>(), GREAT);
    indicator.correctBoundaryConditions();
}


template<class Type>
bool Foam::functionObjects::blendingFactor::calcScheme()
{
    typedef GeometricField<Type, fvPatchField, volMesh> FieldType;

    if (!foundObject<FieldType>(fieldName_, false))
    {
        return false;
    }

    const FieldType& field = lookupObject<FieldType>(fieldName_);

    const word divScheme("div(" + phiName_ + ',' + fieldName_ + ')');
    ITstream& its = mesh_.divScheme(divScheme);

    const surfaceScalarField& phi = lookupObject<surfaceScalarField>(phiName_);

    tmp<fv::convectionScheme<Type>> tcs =
        fv::convectionScheme<Type>::New(mesh_, phi, its);

    if (isA<fv::boundedConvectionScheme<Type>>(tcs()))
    {
        const fv::boundedConvectionScheme<Type>& bcs =
            refCast<const fv::boundedConvectionScheme<Type>>(tcs());

        calcBlendingFactor(field, bcs.scheme());
    }
    else
    {
        const fv::gaussConvectionScheme<Type>& gcs =
            refCast<const fv::gaussConvectionScheme<Type>>(tcs());

        calcBlendingFactor(field, gcs);
    }

    return true;
}


// ************************************************************************* //
