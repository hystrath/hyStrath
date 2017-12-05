/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015-2016 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify i
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

#include "externalCoupled.H"
//#include "fvMesh.H"
#include "OSspecific.H"
#include "IFstream.H"
#include "OFstream.H"
#include "volFields.H"
#include "externalCoupledMixedFvPatchFields.H"
#include "mixedFvPatchFields.H"
#include "fixedGradientFvPatchFields.H"
#include "fixedValueFvPatchFields.H"
#include "OStringStream.H"
#include "globalIndex.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
bool Foam::functionObjects::externalCoupled::readData
(
    const UPtrList<const fvMesh>& meshes,
    const wordRe& groupName,
    const word& fieldName
)
{
    typedef GeometricField<Type, fvPatchField, volMesh> volFieldType;
    typedef externalCoupledMixedFvPatchField<Type> patchFieldType;

    wordList regionNames(meshes.size());
    forAll(meshes, i)
    {
        regionNames[i] = meshes[i].dbDir();
    }

    // File only opened on master; contains data for all processors, for all
    // patchIDs.
    autoPtr<IFstream> masterFilePtr;
    if (Pstream::master())
    {
        const fileName transferFile
        (
            groupDir(commsDir_, compositeName(regionNames), groupName)
          / fieldName + ".in"
        );

        Log << type() << ": reading data from " << transferFile << endl;

        masterFilePtr.reset(new IFstream(transferFile));

        if (!masterFilePtr().good())
        {
            FatalIOErrorInFunction(masterFilePtr())
                << "Cannot open file for region " << compositeName(regionNames)
                << ", field " << fieldName
                << exit(FatalIOError);
        }
    }


    label nFound = 0;
    forAll(meshes, i)
    {
        const fvMesh& mesh = meshes[i];

        if (!mesh.foundObject<volFieldType>(fieldName))
        {
            continue;
        }

        nFound++;

        const volFieldType& cvf = mesh.lookupObject<volFieldType>(fieldName);
        const typename volFieldType::Boundary& bf = cvf.boundaryField();


        // Get the patches
        const labelList patchIDs
        (
            mesh.boundaryMesh().patchSet
            (
                List<wordRe>(1, groupName)
            ).sortedToc()
        );

        // Handle column-wise reading of patch data. Supports most easy types
        forAll(patchIDs, i)
        {
            label patchi = patchIDs[i];

            if (isA<patchFieldType>(bf[patchi]))
            {
                // Explicit handling of externalCoupledMixed bcs - they
                // have specialised reading routines.

                patchFieldType& pf = const_cast<patchFieldType&>
                (
                    refCast<const patchFieldType>
                    (
                        bf[patchi]
                    )
                );

                // Read from master into local stream
                OStringStream os;
                readLines
                (
                    bf[patchi].size(),      // number of lines to read
                    masterFilePtr,
                    os
                );

                // Pass responsability for all reading over to bc
                pf.readData(IStringStream(os.str())());

                // Update the value from the read coefficicient. Bypass any
                // additional processing by derived type.
                pf.patchFieldType::evaluate();
            }
            else if (isA<mixedFvPatchField<Type>>(bf[patchi]))
            {
                // Read columns from file for
                // value, snGrad, refValue, refGrad, valueFraction
                List<scalarField> data;
                readColumns
                (
                    bf[patchi].size(),              // number of lines to read
                    4*pTraits<Type>::nComponents+1, // nColumns: 4*Type+1*scalar
                    masterFilePtr,
                    data
                );

                mixedFvPatchField<Type>& pf =
                const_cast<mixedFvPatchField<Type>&>
                (
                    refCast<const mixedFvPatchField<Type>>
                    (
                        bf[patchi]
                    )
                );

                // Transfer read data to bc.
                // Skip value, snGrad
                direction columni = 2*pTraits<Type>::nComponents;

                Field<Type>& refValue = pf.refValue();
                for
                (
                    direction cmpt = 0;
                    cmpt < pTraits<Type>::nComponents;
                    cmpt++
                )
                {
                    refValue.replace(cmpt, data[columni++]);
                }
                Field<Type>& refGrad = pf.refGrad();
                for
                (
                    direction cmpt = 0;
                    cmpt < pTraits<Type>::nComponents;
                    cmpt++
                )
                {
                    refGrad.replace(cmpt, data[columni++]);
                }
                pf.valueFraction() = data[columni];

                // Update the value from the read coefficicient. Bypass any
                // additional processing by derived type.
                pf.mixedFvPatchField<Type>::evaluate();
            }
            else if (isA<fixedGradientFvPatchField<Type>>(bf[patchi]))
            {
                // Read columns for value and gradient
                List<scalarField> data;
                readColumns
                (
                    bf[patchi].size(),              // number of lines to read
                    2*pTraits<Type>::nComponents,   // nColumns: Type
                    masterFilePtr,
                    data
                );

                fixedGradientFvPatchField<Type>& pf =
                const_cast<fixedGradientFvPatchField<Type>&>
                (
                    refCast<const fixedGradientFvPatchField<Type>>
                    (
                        bf[patchi]
                    )
                );

                // Transfer gradient to bc
                Field<Type>& gradient = pf.gradient();
                for
                (
                    direction cmpt = 0;
                    cmpt < pTraits<Type>::nComponents;
                    cmpt++
                )
                {
                    gradient.replace
                    (
                        cmpt,
                        data[pTraits<Type>::nComponents+cmpt]
                    );
                }

                // Update the value from the read coefficicient. Bypass any
                // additional processing by derived type.
                pf.fixedGradientFvPatchField<Type>::evaluate();
            }
            else if (isA<fixedValueFvPatchField<Type>>(bf[patchi]))
            {
                // Read columns for value only
                List<scalarField> data;
                readColumns
                (
                    bf[patchi].size(),              // number of lines to read
                    pTraits<Type>::nComponents,     // number of columns to read
                    masterFilePtr,
                    data
                );

                // Transfer read value to bc
                Field<Type> value(bf[patchi].size());
                for
                (
                    direction cmpt = 0;
                    cmpt < pTraits<Type>::nComponents;
                    cmpt++
                )
                {
                    value.replace(cmpt, data[cmpt]);
                }

                fixedValueFvPatchField<Type>& pf =
                const_cast<fixedValueFvPatchField<Type>&>
                (
                    refCast<const fixedValueFvPatchField<Type>>
                    (
                        bf[patchi]
                    )
                );

                pf == value;

                // Update the value from the read coefficicient. Bypass any
                // additional processing by derived type.
                pf.fixedValueFvPatchField<Type>::evaluate();
            }
            else
            {
                FatalErrorInFunction
                    << "Unsupported boundary condition " << bf[patchi].type()
                    << " for patch " << bf[patchi].patch().name()
                    << " in region " << mesh.name()
                    << exit(FatalError);
            }

            initialised_ = true;
        }
    }

    return nFound > 0;
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::functionObjects::externalCoupled::gatherAndCombine
(
    const Field<Type>& fld
)
{
    // Collect values from all processors
    List<Field<Type>> gatheredValues(Pstream::nProcs());
    gatheredValues[Pstream::myProcNo()] = fld;
    Pstream::gatherList(gatheredValues);


    tmp<Field<Type>> tresult(new Field<Type>(0));
    Field<Type>& result = tresult.ref();

    if (Pstream::master())
    {
        // Combine values into single field
        label globalElemi = 0;

        forAll(gatheredValues, lsti)
        {
            globalElemi += gatheredValues[lsti].size();
        }

        result.setSize(globalElemi);

        globalElemi = 0;

        forAll(gatheredValues, lsti)
        {
            const Field<Type>& sub = gatheredValues[lsti];

            forAll(sub, elemi)
            {
                result[globalElemi++] = sub[elemi];
            }
        }
    }

    return tresult;
}


template<class Type>
bool Foam::functionObjects::externalCoupled::writeData
(
    const UPtrList<const fvMesh>& meshes,
    const wordRe& groupName,
    const word& fieldName
) const
{
    typedef GeometricField<Type, fvPatchField, volMesh> volFieldType;
    typedef externalCoupledMixedFvPatchField<Type> patchFieldType;

    wordList regionNames(meshes.size());
    forAll(meshes, i)
    {
        regionNames[i] = meshes[i].dbDir();
    }

    // File only opened on master; contains data for all processors, for all
    // patchIDs
    autoPtr<OFstream> masterFilePtr;
    if (Pstream::master())
    {
        const fileName transferFile
        (
            groupDir(commsDir_, compositeName(regionNames), groupName)
          / fieldName + ".out"
        );

        Log << type() << ": writing data to " << transferFile << endl;

        masterFilePtr.reset(new OFstream(transferFile));

        if (!masterFilePtr().good())
        {
            FatalIOErrorInFunction(masterFilePtr())
                << "Cannot open file for region " << compositeName(regionNames)
                << ", field " << fieldName
                << exit(FatalIOError);
        }
    }


    bool headerDone = false;

    label nFound = 0;

    forAll(meshes, i)
    {
        const fvMesh& mesh = meshes[i];

        if (!mesh.foundObject<volFieldType>(fieldName))
        {
            continue;
        }

        nFound++;

        const volFieldType& cvf = mesh.lookupObject<volFieldType>(fieldName);
        const typename volFieldType::Boundary& bf = cvf.boundaryField();


        // Get the patches
        const labelList patchIDs
        (
            mesh.boundaryMesh().patchSet
            (
                List<wordRe>(1, groupName)
            ).sortedToc()
        );

        // Handle column-wise writing of patch data. Supports most easy types
        forAll(patchIDs, i)
        {
            label patchi = patchIDs[i];

            const globalIndex globalFaces(bf[patchi].size());

            if (isA<patchFieldType>(bf[patchi]))
            {
                // Explicit handling of externalCoupledMixed bcs - they
                // have specialised writing routines

                const patchFieldType& pf = refCast<const patchFieldType>
                (
                    bf[patchi]
                );
                OStringStream os;

                // Pass responsibility for all writing over to bc
                pf.writeData(os);

                // Collect contributions from all processors and output them on
                // master
                if (Pstream::master())
                {
                    // Output master data first
                    if (!headerDone)
                    {
                        pf.writeHeader(masterFilePtr());
                        headerDone = true;
                    }
                    masterFilePtr() << os.str().c_str();

                    for (label proci = 1; proci < Pstream::nProcs(); proci++)
                    {
                        IPstream fromSlave(Pstream::scheduled, proci);
                        string str(fromSlave);
                        masterFilePtr() << str.c_str();
                    }
                }
                else
                {
                    OPstream toMaster(Pstream::scheduled, Pstream::masterNo());
                    toMaster << os.str();
                }
            }
            else if (isA<mixedFvPatchField<Type>>(bf[patchi]))
            {
                const mixedFvPatchField<Type>& pf =
                    refCast<const mixedFvPatchField<Type>>(bf[patchi]);

                Field<Type> value(gatherAndCombine(pf));
                Field<Type> snGrad(gatherAndCombine(pf.snGrad()()));
                Field<Type> refValue(gatherAndCombine(pf.refValue()));
                Field<Type> refGrad(gatherAndCombine(pf.refGrad()));
                scalarField valueFraction(gatherAndCombine(pf.valueFraction()));

                if (Pstream::master())
                {
                    forAll(refValue, facei)
                    {
                        masterFilePtr()
                            << value[facei] << token::SPACE
                            << snGrad[facei] << token::SPACE
                            << refValue[facei] << token::SPACE
                            << refGrad[facei] << token::SPACE
                            << valueFraction[facei] << nl;
                    }
                }
            }
            else
            {
                // Output the value and snGrad
                Field<Type> value(gatherAndCombine(bf[patchi]));
                Field<Type> snGrad(gatherAndCombine(bf[patchi].snGrad()()));
                if (Pstream::master())
                {
                    forAll(value, facei)
                    {
                        masterFilePtr()
                            << value[facei] << token::SPACE
                            << snGrad[facei] << nl;
                    }
                }
            }
        }
    }

    return nFound > 0;
}


// ************************************************************************* //
