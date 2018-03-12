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

#include "mapFields.H"
#include "meshToMesh.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(mapFields, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        mapFields,
        dictionary
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::functionObjects::mapFields::createInterpolation
(
    const dictionary& dict
)
{
    const fvMesh& meshTarget = mesh_;
    const word mapRegionName(dict.lookup("mapRegion"));

    Info<< name() << ':' << nl
        << "    Reading mesh " << mapRegionName << endl;

    mapRegionPtr_.reset
    (
        new fvMesh
        (
            IOobject
            (
                mapRegionName,
                meshTarget.time().constant(),
                meshTarget.time()
            )
        )
    );
    const fvMesh& mapRegion = mapRegionPtr_();
    word mapMethodName(dict.lookup("mapMethod"));
    if (!meshToMesh::interpolationMethodNames_.found(mapMethodName))
    {
        FatalErrorInFunction
            << type() << " " << name() << ": unknown map method "
            << mapMethodName << nl
            << "Available methods include: "
            << meshToMesh::interpolationMethodNames_.sortedToc()
            << exit(FatalError);
    }

    meshToMesh::interpolationMethod mapMethod
    (
        meshToMesh::interpolationMethodNames_[mapMethodName]
    );

    // Lookup corresponding AMi method
    word patchMapMethodName =
        AMIPatchToPatchInterpolation::interpolationMethodToWord
        (
            meshToMesh::interpolationMethodAMI(mapMethod)
        );

    // Optionally override
    if (dict.readIfPresent("patchMapMethod", patchMapMethodName))
    {
        Info<< "    Patch mapping method: " << patchMapMethodName << endl;
    }

    bool consistent = readBool(dict.lookup("consistent"));

    Info<< "    Creating mesh to mesh interpolation" << endl;

    if (consistent)
    {
        interpPtr_.reset
        (
            new meshToMesh
            (
                mapRegion,
                meshTarget,
                mapMethodName,
                patchMapMethodName
            )
        );
    }
    else
    {
        HashTable<word> patchMap(dict.lookup("patchMap"));
        wordList cuttingPatches(dict.lookup("cuttingPatches"));

        interpPtr_.reset
        (
            new meshToMesh
            (
                mapRegion,
                meshTarget,
                mapMethodName,
                patchMapMethodName,
                patchMap,
                cuttingPatches
            )
        );
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::mapFields::mapFields
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    mapRegionPtr_(),
    interpPtr_(),
    fieldNames_()
{
    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::mapFields::~mapFields()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::mapFields::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);

    dict.lookup("fields") >> fieldNames_;
    createInterpolation(dict);
    return true;
}


bool Foam::functionObjects::mapFields::execute()
{
    return true;
}


bool Foam::functionObjects::mapFields::write()
{
    Log << type() << " " << name() << " write:" << nl;

    bool ok = false;

    ok = writeFieldType<scalar>() || ok;
    ok = writeFieldType<vector>() || ok;
    ok = writeFieldType<sphericalTensor>() || ok;
    ok = writeFieldType<symmTensor>() || ok;
    ok = writeFieldType<tensor>() || ok;

    if (log)
    {
        if (!ok)
        {
            Info<< "    none" << nl;
        }

        Info<< endl;
    }

    return true;
}


// ************************************************************************* //
