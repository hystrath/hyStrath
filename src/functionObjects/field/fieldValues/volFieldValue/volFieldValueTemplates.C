/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2015-2016 OpenCFD Ltd.
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

#include "volFieldValue.H"
#include "volFields.H"

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class Type>
bool Foam::functionObjects::fieldValues::volFieldValue::validField
(
    const word& fieldName
) const
{
    typedef GeometricField<Type, fvPatchField, volMesh> vf;

    if (obr_.foundObject<vf>(fieldName))
    {
        return true;
    }

    return false;
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::functionObjects::fieldValues::volFieldValue::getFieldValues
(
    const word& fieldName,
    const bool mustGet
) const
{
    typedef GeometricField<Type, fvPatchField, volMesh> vf;

    if (obr_.foundObject<vf>(fieldName))
    {
        return filterField(obr_.lookupObject<vf>(fieldName));
    }

    if (mustGet)
    {
        FatalErrorInFunction
            << "Field " << fieldName << " not found in database"
            << abort(FatalError);
    }

    return tmp<Field<Type>>(new Field<Type>(0.0));
}


template<class Type>
Type Foam::functionObjects::fieldValues::volFieldValue::processValues
(
    const Field<Type>& values,
    const scalarField& V,
    const scalarField& weightField
) const
{
    Type result = Zero;
    switch (operation_)
    {
        case opSum:
        {
            result = gSum(values);
            break;
        }
        case opSumMag:
        {
            result = gSum(cmptMag(values));
            break;
        }
        case opAverage:
        {
            label n = returnReduce(values.size(), sumOp<label>());
            result = gSum(values)/(scalar(n) + ROOTVSMALL);
            break;
        }
        case opWeightedAverage:
        {
            label wSize = returnReduce(weightField.size(), sumOp<label>());

            if (wSize > 0)
            {
                result =
                    gSum(weightField*values)/(gSum(weightField) + ROOTVSMALL);
            }
            else
            {
                label n = returnReduce(values.size(), sumOp<label>());
                result = gSum(values)/(scalar(n) + ROOTVSMALL);
            }
            break;
        }
        case opVolAverage:
        {
            result = gSum(values*V)/(gSum(V) + ROOTVSMALL);
            break;
        }
        case opWeightedVolAverage:
        {
            result = gSum(weightField*V*values)/gSum(weightField*V);
            break;
        }
        case opVolIntegrate:
        {
            result = gSum(V*values);
            break;
        }
        case opMin:
        {
            result = gMin(values);
            break;
        }
        case opMax:
        {
            result = gMax(values);
            break;
        }
        case opCoV:
        {
            const scalar sumV = gSum(V);

            Type meanValue = gSum(V*values)/sumV;

            const label nComp = pTraits<Type>::nComponents;

            for (direction d=0; d<nComp; ++d)
            {
                scalarField vals(values.component(d));
                scalar mean = component(meanValue, d);
                scalar& res = setComponent(result, d);

                res = sqrt(gSum(V*sqr(vals - mean))/sumV)/(mean + ROOTVSMALL);
            }

            break;
        }
        case opNone:
        {}
    }

    return result;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
bool Foam::functionObjects::fieldValues::volFieldValue::writeValues
(
    const word& fieldName,
    const scalarField& weightField
)
{
    const bool ok = validField<Type>(fieldName);

    if (ok)
    {
        Field<Type> values(getFieldValues<Type>(fieldName));
        scalarField V(filterField(fieldValue::mesh_.V()));

        if (writeFields_)
        {
            Field<Type> allValues(values);
            combineFields(allValues);

            if (Pstream::master())
            {
                IOField<Type>
                (
                    IOobject
                    (
                        fieldName + '_' + regionTypeNames_[regionType_]
                      + '-' + volRegion::regionName_,
                        obr_.time().timeName(),
                        obr_,
                        IOobject::NO_READ,
                        IOobject::NO_WRITE
                    ),
                    scaleFactor_*weightField*allValues
                ).write();
            }
        }

        // Apply scale factor
        values *= scaleFactor_;

        Type result = processValues(values, V, weightField);

        file()<< tab << result;

        Log << "    " << operationTypeNames_[operation_]
            << "(" << volRegion::regionName_ << ") of " << fieldName
            <<  " = " << result << endl;

        // Write state/results information
        const word& opName = operationTypeNames_[operation_];
        word resultName =
            opName + '(' + volRegion::regionName_ + ',' + fieldName + ')';
        this->setResult(resultName, result);
    }

    return ok;
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::functionObjects::fieldValues::volFieldValue::filterField
(
    const Field<Type>& field
) const
{
    if (isNull(cellIDs()))
    {
        return field;
    }
    else
    {
        return tmp<Field<Type>>(new Field<Type>(field, cellIDs()));
    }
}


// ************************************************************************* //
