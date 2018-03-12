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

#include "fvMesh.H"
#include "polyPatch.H"
#include "lduSchedule.H"
#include "meshToMesh.H"

template<class Type>
void Foam::functionObjects::mapFields::evaluateConstraintTypes
(
    GeometricField<Type, fvPatchField, volMesh>& fld
) const
{
    typedef GeometricField<Type, fvPatchField, volMesh> VolFieldType;

    typename VolFieldType::Boundary& fldBf = fld.boundaryFieldRef();

    if
    (
        Pstream::defaultCommsType == Pstream::blocking
     || Pstream::defaultCommsType == Pstream::nonBlocking
    )
    {
        label nReq = Pstream::nRequests();

        forAll(fldBf, patchi)
        {
            fvPatchField<Type>& tgtField = fldBf[patchi];

            if
            (
                tgtField.type() == tgtField.patch().patch().type()
             && polyPatch::constraintType(tgtField.patch().patch().type())
            )
            {
                tgtField.initEvaluate(Pstream::defaultCommsType);
            }
        }

        // Block for any outstanding requests
        if
        (
            Pstream::parRun()
         && Pstream::defaultCommsType == Pstream::nonBlocking
        )
        {
            Pstream::waitRequests(nReq);
        }

        forAll(fldBf, patchi)
        {
            fvPatchField<Type>& tgtField = fldBf[patchi];

            if
            (
                tgtField.type() == tgtField.patch().patch().type()
             && polyPatch::constraintType(tgtField.patch().patch().type())
            )
            {
                tgtField.evaluate(Pstream::defaultCommsType);
            }
        }
    }
    else if (Pstream::defaultCommsType == Pstream::scheduled)
    {
        const lduSchedule& patchSchedule =
            fld.mesh().globalData().patchSchedule();

        forAll(patchSchedule, patchEvali)
        {
            label patchi = patchSchedule[patchEvali].patch;
            fvPatchField<Type>& tgtField = fldBf[patchi];

            if
            (
                tgtField.type() == tgtField.patch().patch().type()
             && polyPatch::constraintType(tgtField.patch().patch().type())
            )
            {
                if (patchSchedule[patchEvali].init)
                {
                    tgtField.initEvaluate(Pstream::scheduled);
                }
                else
                {
                    tgtField.evaluate(Pstream::scheduled);
                }
            }
        }
    }
}


template<class Type>
bool Foam::functionObjects::mapFields::writeFieldType() const
{
    typedef GeometricField<Type, fvPatchField, volMesh> VolFieldType;

    const fvMesh& mapRegion = mapRegionPtr_();

    wordList fieldNames(this->mesh_.names(VolFieldType::typeName));
    labelList selected = findStrings(fieldNames_, fieldNames);
    forAll(selected, i)
    {
        const word& fieldName = fieldNames[selected[i]];
        const VolFieldType& field = lookupObject<VolFieldType>(fieldName);

        Log << "    " << fieldName;

        IOobject mapRegionIO
        (
            fieldName,
            time_.timeName(),
            mapRegion
        );

        tmp<VolFieldType> tfieldMapRegion(interpPtr_->mapTgtToSrc(field));

        Log << ": interpolated";

        VolFieldType fieldMapRegion(mapRegionIO, tfieldMapRegion);

        evaluateConstraintTypes(fieldMapRegion);

        fieldMapRegion.write();

        Log << " and written" << nl;
    }

    return selected.size() > 0;
}


// ************************************************************************* //
