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

Class
    Foam::functionObjects::mapFields

Group
    grpFieldFunctionObjects

Description
    Map fields from local mesh to secondary mesh at run-time.

Usage
    Example of function object specification to map fields:
    \verbatim
    mapFields1
    {
        type            mapFields;
        libs            ("libfieldFunctionObjects.so");
        ...

        mapRegion       coarseMesh;
        mapMethod       cellVolumeWeight;
        consistent      yes;

        fields          ("U.*" p);
    }

    Where the entries comprise:
    \table
        Property     | Description                | Required  | Default value
        type         | Type name: mapFields       | yes       |
        mapRegion    | Name of region to map to   | yes       |
        mapMethod    | Mapping method             | yes       |
        patchMapMethod | Patch mapping method     | no        | <auto>
        consistent   | Mapping meshes have consistent boundaries | yes |
        fields       | List of field names to map | yes       |
        log          | Log to standard output     | no        | yes
    \endtable


SourceFiles
    mapFields.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_mapFields_H
#define functionObjects_mapFields_H

#include "fvMeshFunctionObject.H"
#include "volFieldsFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class meshToMesh;

namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                          Class mapFields Declaration
\*---------------------------------------------------------------------------*/

class mapFields
:
    public fvMeshFunctionObject
{
    // Private data

        //- Locally cached map region mesh (map to this mesh)
        autoPtr<fvMesh> mapRegionPtr_;

        //- Mesh-to-mesh interpolation
        autoPtr<meshToMesh> interpPtr_;

        //- List of field names to interpolate
        wordReList fieldNames_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        mapFields(const mapFields&) = delete;

        //- Disallow default bitwise assignment
        void operator=(const mapFields&) = delete;

        //- Helper function to create the mesh-to-mesh interpolation
        void createInterpolation(const dictionary& dict);

        //- Helper function to evaluate constraint patches after mapping
        template<class Type>
        void evaluateConstraintTypes
        (
            GeometricField<Type, fvPatchField, volMesh>& fld
        ) const;

        //- Helper function to interpolate and write the field
        template<class Type>
        bool writeFieldType() const;


public:

    //- Runtime type information
    TypeName("mapFields");


    // Constructors

        //- Construct from Time and dictionary
        mapFields
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );


    //- Destructor
    virtual ~mapFields();


    // Member Functions

        //- Read the mapFields data
        virtual bool read(const dictionary&);

        //- Execute, currently does nothing
        virtual bool execute();

        //- Calculate the mapFields and write
        virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "mapFieldsTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
