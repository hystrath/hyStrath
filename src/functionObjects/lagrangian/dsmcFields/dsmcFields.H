/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
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
    Foam::functionObjects::dsmcFields

Group
    grpLagrangianFunctionObjects

Description
    This function object calculates and outputs the intensive fields:
    - UMean
    - translationalT
    - internalT
    - overallT
    from averaged extensive fields from a DSMC calculation.

Usage
    Example of function object specification to calculate DSMC fields:
    \verbatim
    dsmcFields1
    {
        type        dsmcFields;
        libs        ("liblagrangianFunctionObjects.so");
        ...
    }
    \endverbatim

    Where the entries comprise:
    \table
        Property     | Description             | Required    | Default value
        type         | Type name: dsmcFields   | yes         |
        log          | Log to standard output  | no          | yes
    \endtable


SourceFiles
    dsmcFields.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_dsmcFields_H
#define functionObjects_dsmcFields_H

#include "fvMeshFunctionObject.H"
#include "Switch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                         Class dsmcFields Declaration
\*---------------------------------------------------------------------------*/

class dsmcFields
:
    public fvMeshFunctionObject
{

        //- Switch to send output to Info as well as to file
        Switch log_;
    // Private Member Functions

        //- Disallow default bitwise copy construct
        dsmcFields(const dsmcFields&);

        //- Disallow default bitwise assignment
        void operator=(const dsmcFields&);


public:

    //- Runtime type information
    TypeName("dsmcFields");


    // Constructors

        //- Construct from Time and dictionary
        dsmcFields
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );


    //- Destructor
    virtual ~dsmcFields();


    // Member Functions

        //- Read the dsmcFields data
        virtual bool read(const dictionary&);

        //- Do nothing
        virtual bool execute();

        //- Calculate and write the DSMC fields
        virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
