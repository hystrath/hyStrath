/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013-2015 OpenFOAM Foundation
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

Class
    Foam::functionObjects::DESModelRegions

Group
    grpFieldFunctionObjects

Description
    This function object writes out an indicator field for DES turbulence
    calculations, that is:
    - 0 for RAS regions
    - 1 for LES regions

    The field is stored on the mesh database so that it can be retrieved and
    used for other applications.

Usage
    Example of function object specification to generate DES indicator field:
    \verbatim
    DESModelRegions1
    {
        type        DESModelRegions;
        libs        ("libfieldFunctionObjects.so");
        ...
    }
    \endverbatim

    Where the entries comprise:
    \table
        Property     | Description             | Required    | Default value
        type         | type name: DESModelRegions| yes       |
        resultName   | Name of DES indicator field | no      | <function name>
        log          | log to standard output  | no          | yes
    \endtable

SourceFiles
    DESModelRegions.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_DESModelRegions_H
#define functionObjects_DESModelRegions_H

#include "fvMeshFunctionObject.H"
#include "writeFile.H"
#include "volFieldsFwd.H"
#include "OFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                       Class DESModelRegions Declaration
\*---------------------------------------------------------------------------*/

class DESModelRegions
:
    public fvMeshFunctionObject,
    public writeFile
{
protected:

    // Protected data

        //- Result name
        word resultName_;


    // Protected Member Functions

        //- File header information
        virtual void writeFileHeader(Ostream& os) const;

        //- Disallow default bitwise copy construct
        DESModelRegions(const DESModelRegions&) = delete;

        //- Disallow default bitwise assignment
        void operator=(const DESModelRegions&) = delete;


public:

    //- Runtime type information
    TypeName("DESModelRegions");


    // Constructors
        //- Construct from Time and dictionary
        DESModelRegions
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );


    //- Destructor
    virtual ~DESModelRegions();


    // Member Functions

        //- Read the DESModelRegions data
        virtual bool read(const dictionary&);

        //- Execute
        virtual bool execute();

        //- Calculate the DESModelRegions and write
        virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
