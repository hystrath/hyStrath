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

#include "fieldAverageItem.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "surfFields.H"
#include "OFstream.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
void Foam::functionObjects::fieldAverage::addMeanFieldType(const label fieldi)
{
    // Field has been found, so set active flag to true
    faItems_[fieldi].active() = true;

    const word& fieldName = faItems_[fieldi].fieldName();
    const word& meanFieldName = faItems_[fieldi].meanFieldName();

    Log << "    Reading/initialising field " << meanFieldName << endl;

    if (foundObject<Type>(meanFieldName))
    {
       // do nothing
    }
    else if (obr().found(meanFieldName))
    {
        Log << "    Cannot allocate average field " << meanFieldName
            << " since an object with that name already exists."
            << " Disabling averaging for field." << endl;

        faItems_[fieldi].mean() = false;
    }
    else
    {
        const Type& baseField = lookupObject<Type>(fieldName);

        // Store on registry
        obr().store
        (
            new Type
            (
                IOobject
                (
                    meanFieldName,
                    obr().time().timeName(obr().time().startTime().value()),
                    obr(),
                    restartOnOutput_
                  ? IOobject::NO_READ
                  : IOobject::READ_IF_PRESENT,
                    IOobject::NO_WRITE
                ),
                1*baseField
            )
        );
    }
}


template<class Type>
void Foam::functionObjects::fieldAverage::addMeanField(const label fieldi)
{
    typedef GeometricField<Type, fvPatchField, volMesh> VolFieldType;
    typedef GeometricField<Type, fvsPatchField, surfaceMesh> SurfaceFieldType;
    typedef DimensionedField<Type, surfGeoMesh> SurfFieldType;

    if (faItems_[fieldi].mean())
    {
        const word& fieldName = faItems_[fieldi].fieldName();

        if (foundObject<VolFieldType>(fieldName))
        {
            addMeanFieldType<VolFieldType>(fieldi);
        }
        else if (foundObject<SurfaceFieldType>(fieldName))
        {
            addMeanFieldType<SurfaceFieldType>(fieldi);
        }
        else if (foundObject<SurfFieldType>(fieldName))
        {
            addMeanFieldType<SurfFieldType>(fieldi);
        }
    }
}


template<class Type1, class Type2>
void Foam::functionObjects::fieldAverage::addPrime2MeanFieldType
(
    const label fieldi
)
{
    const word& fieldName = faItems_[fieldi].fieldName();
    const word& meanFieldName = faItems_[fieldi].meanFieldName();
    const word& prime2MeanFieldName = faItems_[fieldi].prime2MeanFieldName();

    Log << "    Reading/initialising field " << prime2MeanFieldName << nl;

    if (foundObject<Type2>(prime2MeanFieldName))
    {
        // do nothing
    }
    else if (obr().found(prime2MeanFieldName))
    {
        Log << "    Cannot allocate average field " << prime2MeanFieldName
            << " since an object with that name already exists."
            << " Disabling averaging for field." << endl;

        faItems_[fieldi].prime2Mean() = false;
    }
    else
    {
        const Type1& baseField = lookupObject<Type1>(fieldName);
        const Type1& meanField = lookupObject<Type1>(meanFieldName);

        // Store on registry
        obr().store
        (
            new Type2
            (
                IOobject
                (
                    prime2MeanFieldName,
                    obr().time().timeName(obr().time().startTime().value()),
                    obr(),
                    restartOnOutput_
                  ? IOobject::NO_READ
                  : IOobject::READ_IF_PRESENT,
                    IOobject::NO_WRITE
                ),
                sqr(baseField) - sqr(meanField)
            )
        );
    }
}


template<class Type1, class Type2>
void Foam::functionObjects::fieldAverage::addPrime2MeanField(const label fieldi)
{
    typedef GeometricField<Type1, fvPatchField, volMesh> VolFieldType1;
    typedef GeometricField<Type1, fvsPatchField, surfaceMesh> SurfaceFieldType1;
    typedef DimensionedField<Type1, surfGeoMesh> SurfFieldType1;

    typedef GeometricField<Type2, fvPatchField, volMesh> VolFieldType2;
    typedef GeometricField<Type2, fvsPatchField, surfaceMesh> SurfaceFieldType2;
    typedef DimensionedField<Type2, surfGeoMesh> SurfFieldType2;

    if (faItems_[fieldi].prime2Mean())
    {
        const word& fieldName = faItems_[fieldi].fieldName();

        if (!faItems_[fieldi].mean())
        {
            FatalErrorInFunction
                << "To calculate the prime-squared average, the "
                << "mean average must also be selected for field "
                << fieldName << nl << exit(FatalError);
        }

        if (foundObject<VolFieldType1>(fieldName))
        {
            addPrime2MeanFieldType<VolFieldType1, VolFieldType2>(fieldi);
        }
        else if (foundObject<SurfaceFieldType1>(fieldName))
        {
            addPrime2MeanFieldType<SurfaceFieldType1, SurfaceFieldType2>
            (
                fieldi
            );
        }
        else if (foundObject<SurfFieldType1>(fieldName))
        {
            addPrime2MeanFieldType<SurfFieldType1, SurfFieldType2>
            (
                fieldi
            );
        }
    }
}


template<class Type>
void Foam::functionObjects::fieldAverage::calculateMeanFieldType
(
    const label fieldi
) const
{
    const word& fieldName = faItems_[fieldi].fieldName();

    if (foundObject<Type>(fieldName))
    {
        const Type& baseField = lookupObject<Type>(fieldName);

        Type& meanField = const_cast<Type&>
        (
            lookupObject<Type>(faItems_[fieldi].meanFieldName())
        );

        scalar dt = obr().time().deltaTValue();
        scalar Dt = totalTime_[fieldi];

        if (faItems_[fieldi].iterBase())
        {
            dt = 1;
            Dt = scalar(totalIter_[fieldi]);
        }

        scalar beta = dt/Dt;

        if (faItems_[fieldi].window() > 0)
        {
            const scalar w = faItems_[fieldi].window();

            if (Dt - dt >= w)
            {
                beta = dt/w;
            }
        }

        meanField = (1 - beta)*meanField + beta*baseField;
    }
}


template<class Type>
void Foam::functionObjects::fieldAverage::calculateMeanFields() const
{
    typedef GeometricField<Type, fvPatchField, volMesh> VolFieldType;
    typedef GeometricField<Type, fvsPatchField, surfaceMesh> SurfaceFieldType;
    typedef DimensionedField<Type, surfGeoMesh> SurfFieldType;

    forAll(faItems_, i)
    {
        if (faItems_[i].mean())
        {
            calculateMeanFieldType<VolFieldType>(i);
            calculateMeanFieldType<SurfaceFieldType>(i);
            calculateMeanFieldType<SurfFieldType>(i);
        }
    }
}


template<class Type1, class Type2>
void Foam::functionObjects::fieldAverage::calculatePrime2MeanFieldType
(
    const label fieldi
) const
{
    const word& fieldName = faItems_[fieldi].fieldName();

    if (foundObject<Type1>(fieldName))
    {
        const Type1& baseField = lookupObject<Type1>(fieldName);
        const Type1& meanField =
            lookupObject<Type1>(faItems_[fieldi].meanFieldName());

        Type2& prime2MeanField =
            lookupObjectRef<Type2>(faItems_[fieldi].prime2MeanFieldName());

        scalar dt = obr().time().deltaTValue();
        scalar Dt = totalTime_[fieldi];

        if (faItems_[fieldi].iterBase())
        {
            dt = 1;
            Dt = scalar(totalIter_[fieldi]);
        }

        scalar beta = dt/Dt;

        if (faItems_[fieldi].window() > 0)
        {
            const scalar w = faItems_[fieldi].window();

            if (Dt - dt >= w)
            {
                beta = dt/w;
            }
        }

        prime2MeanField =
            (1 - beta)*prime2MeanField
          + beta*sqr(baseField)
          - sqr(meanField);
    }
}


template<class Type1, class Type2>
void Foam::functionObjects::fieldAverage::calculatePrime2MeanFields() const
{
    typedef GeometricField<Type1, fvPatchField, volMesh> VolFieldType1;
    typedef GeometricField<Type1, fvsPatchField, surfaceMesh> SurfaceFieldType1;
    typedef DimensionedField<Type1, surfGeoMesh> SurfFieldType1;

    typedef GeometricField<Type2, fvPatchField, volMesh> VolFieldType2;
    typedef GeometricField<Type2, fvsPatchField, surfaceMesh> SurfaceFieldType2;
    typedef DimensionedField<Type2, surfGeoMesh> SurfFieldType2;

    forAll(faItems_, i)
    {
        if (faItems_[i].prime2Mean())
        {
            calculatePrime2MeanFieldType<VolFieldType1, VolFieldType2>(i);
            calculatePrime2MeanFieldType<SurfaceFieldType1, SurfaceFieldType2>
            (
                i
            );

            calculatePrime2MeanFieldType<SurfFieldType1, SurfFieldType2>
            (
                i
            );
        }
    }
}


template<class Type1, class Type2>
void Foam::functionObjects::fieldAverage::addMeanSqrToPrime2MeanType
(
    const label fieldi
) const
{
    const word& fieldName = faItems_[fieldi].fieldName();

    if (foundObject<Type1>(fieldName))
    {
        const Type1& meanField =
            lookupObject<Type1>(faItems_[fieldi].meanFieldName());

        Type2& prime2MeanField =
            lookupObjectRef<Type2>(faItems_[fieldi].prime2MeanFieldName());

        prime2MeanField += sqr(meanField);
    }
}


template<class Type1, class Type2>
void Foam::functionObjects::fieldAverage::addMeanSqrToPrime2Mean() const
{
    typedef GeometricField<Type1, fvPatchField, volMesh> VolFieldType1;
    typedef GeometricField<Type1, fvsPatchField, surfaceMesh> SurfaceFieldType1;
    typedef DimensionedField<Type1, surfGeoMesh> SurfFieldType1;

    typedef GeometricField<Type2, fvPatchField, volMesh> VolFieldType2;
    typedef GeometricField<Type2, fvsPatchField, surfaceMesh> SurfaceFieldType2;
    typedef DimensionedField<Type2, surfGeoMesh> SurfFieldType2;

    forAll(faItems_, i)
    {
        if (faItems_[i].prime2Mean())
        {
            addMeanSqrToPrime2MeanType<VolFieldType1, VolFieldType2>(i);
            addMeanSqrToPrime2MeanType<SurfaceFieldType1, SurfaceFieldType2>(i);
            addMeanSqrToPrime2MeanType<SurfFieldType1, SurfFieldType2>(i);
        }
    }
}


template<class Type>
void Foam::functionObjects::fieldAverage::writeFieldType
(
    const word& fieldName
) const
{
    if (foundObject<Type>(fieldName))
    {
        const Type& f = lookupObject<Type>(fieldName);
        f.write();
    }
}


template<class Type>
void Foam::functionObjects::fieldAverage::writeFields() const
{
    typedef GeometricField<Type, fvPatchField, volMesh> VolFieldType;
    typedef GeometricField<Type, fvsPatchField, surfaceMesh> SurfaceFieldType;
    typedef DimensionedField<Type, surfGeoMesh> SurfFieldType;

    forAll(faItems_, i)
    {
        if (faItems_[i].mean())
        {
            const word& fieldName = faItems_[i].meanFieldName();
            writeFieldType<VolFieldType>(fieldName);
            writeFieldType<SurfaceFieldType>(fieldName);
            writeFieldType<SurfFieldType>(fieldName);
        }
        if (faItems_[i].prime2Mean())
        {
            const word& fieldName = faItems_[i].prime2MeanFieldName();
            writeFieldType<VolFieldType>(fieldName);
            writeFieldType<SurfaceFieldType>(fieldName);
            writeFieldType<SurfFieldType>(fieldName);
        }
    }
}


// ************************************************************************* //
