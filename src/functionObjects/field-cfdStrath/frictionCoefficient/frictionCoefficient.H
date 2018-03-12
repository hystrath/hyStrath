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
    Foam::functionObjects::frictionCoefficient

Group
    grpFieldFunctionObjects

Description
    Evaluates and outputs turbulence y+ for turbulence models.

    The field is stored on the mesh database so that it can be
    retrieved and used for other applications.

Usage
    Example of function object specification to calculate the friction
    coefficient:
    \verbatim
    frictionCoefficient1
    {
        type        frictionCoefficient;
        libs        ("libstrathFieldFunctionObjects.so");
        ...
    }
    \endverbatim

    Where the entries comprise:
    \table
        Property     | Description             | Required    | Default value
        type         | Type name: frictionCoefficient  | yes |
        wallShearStress | Name of wallShearStress field | no | wallShearStress
        inflowPatchName | Name of inflow patch | no          | inlet
        result       | Name of frictionCoefficient field | no | \<function name\>
        log          | Log to standard output  | no          | yes
    \endtable


See also
    Foam::functionObject
    Foam::functionObjects::fvMeshFunctionObject
    Foam::functionObjects::writeFile
    Foam::functionObjects::timeControl

SourceFiles
    frictionCoefficient.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_frictionCoefficient_H
#define functionObjects_frictionCoefficient_H

#include "fvMeshFunctionObject.H"
#include "writeFile.H"
#include "volFieldsFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                    Class frictionCoefficient Declaration
\*---------------------------------------------------------------------------*/

class frictionCoefficient
:
    public fvMeshFunctionObject,
    public writeFile
{
    // Private data

        //- Name of wallShearStress field, default is "wallShearStress"
        word wallShearStress_;
        
        //- Name of inflow patch
        word inflowPatchName_;
        
        //- header of wallShearStress field
        IOobject wallShearStressHeader_;
        
    
    // Private Member Functions

        //- File header information
        virtual void writeFileHeader(Ostream& os) const;

        //- Disallow default bitwise copy construct
        frictionCoefficient(const frictionCoefficient&);

        //- Disallow default bitwise assignment
        void operator=(const frictionCoefficient&);


public:

    //- Runtime type information
    TypeName("frictionCoefficient");


    // Constructors

        //- Construct from Time and dictionary
        frictionCoefficient
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );


    //- Destructor
    virtual ~frictionCoefficient();


    // Member Functions

        //- Read the frictionCoefficient data
        virtual bool read(const dictionary&);

        //- Calculate the frictionCoefficient field
        virtual bool execute();

        //- Write the frictionCoefficient field
        virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
