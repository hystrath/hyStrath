/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013-2016 OpenFOAM Foundation
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

Class
    Foam::functionObjects::pressureCoefficient

Group
    grpFieldFunctionObjects

Description
    Evaluates and outputs turbulence y+ for turbulence models.

    The field is stored on the mesh database so that it can be
    retrieved and used for other applications.

Usage
    Example of function object specification to calculate the pressure
    coefficient:
    \verbatim
    pressureCoefficient1
    {
        type        pressureCoefficient;
        libs        ("libstrathFieldFunctionObjects.so");
        ...
    }
    \endverbatim

    Where the entries comprise:
    \table
        Property     | Description             | Required    | Default value
        type         | Type name: pressureCoefficient  | yes |
        inflowPatchName | Name of inflow patch | no          | inlet
        result       | Name of pressureCoefficient field | no | \<function name\>
        log          | Log to standard output  | no          | yes
    \endtable


See also
    Foam::functionObject
    Foam::functionObjects::fvMeshFunctionObject
    Foam::functionObjects::writeFile
    Foam::functionObjects::timeControl

SourceFiles
    pressureCoefficient.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_pressureCoefficient_H
#define functionObjects_pressureCoefficient_H

#include "fvMeshFunctionObject.H"
#include "writeFile.H"
#include "volFieldsFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                    Class pressureCoefficient Declaration
\*---------------------------------------------------------------------------*/

class pressureCoefficient
:
    public fvMeshFunctionObject,
    public writeFile
{
    // Private data

        //- Name of inflow patch
        word inflowPatchName_;
        
    
    // Private Member Functions

        //- File header information
        virtual void writeFileHeader(Ostream& os) const;

        //- Disallow default bitwise copy construct
        pressureCoefficient(const pressureCoefficient&);

        //- Disallow default bitwise assignment
        void operator=(const pressureCoefficient&);


public:

    //- Runtime type information
    TypeName("pressureCoefficient");


    // Constructors

        //- Construct from Time and dictionary
        pressureCoefficient
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );


    //- Destructor
    virtual ~pressureCoefficient();


    // Member Functions

        //- Read the pressureCoefficient data
        virtual bool read(const dictionary&);

        //- Calculate the pressureCoefficient field
        virtual bool execute();

        //- Write the pressureCoefficient field
        virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
