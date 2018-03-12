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

#include "polyPatch.H"
#include "Time.H"
#include "zeroGradientFvPatchField.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
bool Foam::functionObjects::zeroGradient::accept
(
    const GeometricField<Type, fvPatchField, volMesh>& input
)
{
    const typename GeometricField<Type, fvPatchField, volMesh>::Boundary&
        patches = input.boundaryField();

    forAll(patches, patchi)
    {
        if (!polyPatch::constraintType(patches[patchi].patch().patch().type()))
        {
            return true;
        }
    }

    return false;
}


template<class Type>
int Foam::functionObjects::zeroGradient::apply
(
    const word& inputName,
    int& state
)
{
    typedef GeometricField<Type, fvPatchField, volMesh> VolFieldType;

    // State: return 0 (not-processed), -1 (skip), +1 ok

    // Already done, or not available
    if (state || !foundObject<VolFieldType>(inputName))
    {
        return state;
    }

    const VolFieldType& input = lookupObject<VolFieldType>(inputName);

    if (!returnReduce(accept(input), orOp<bool>()))
    {
        state = -1;
        return state;
    }

    word outputName(resultName_);
    outputName.replace("@@", inputName);

    // Also save the field-type, just in case we want it later
    results_.set(outputName, VolFieldType::typeName);

    if (!foundObject<VolFieldType>(outputName))
    {
        tmp<VolFieldType> tzg
        (
            new VolFieldType
            (
                IOobject
                (
                    outputName,
                    time_.timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh_,
                dimensioned<Type>("0", input.dimensions(), Zero),
                zeroGradientFvPatchField<Type>::typeName
            )
        );

        store(outputName, tzg);
    }

    VolFieldType& output =
        const_cast<VolFieldType&>(lookupObject<VolFieldType>(outputName));

    output = input;
    output.correctBoundaryConditions();

    state = +1;
    return state;
}


// ************************************************************************* //
