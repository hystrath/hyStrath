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
    Foam::functionObjects::StantonNo

Group
    grpFieldFunctionObjects

Description
    Evaluates and outputs turbulence y+ for turbulence models.

    The field is stored on the mesh database so that it can be
    retrieved and used for other applications.

Usage
    Example of function object specification to calculate the StantonNo:
    \verbatim
    StantonNo1
    {
        type        StantonNo;
        libs        ("libstrathFieldFunctionObjects.so");
        ...
    }
    \endverbatim

    Where the entries comprise:
    \table
        Property     | Description             | Required    | Default value
        type         | Type name: StantonNo    | yes         |
        wallHeatFlux | Name of wallHeatFlux field | no       | wallHeatFlux
        inflowPatchName | Name of inflow patch | no          | inlet
        result       | Name of Stanton field   | no          | \<function name\>
        log          | Log to standard output  | no          | yes
    \endtable


See also
    Foam::functionObject
    Foam::functionObjects::fvMeshFunctionObject
    Foam::functionObjects::writeFile
    Foam::functionObjects::timeControl

SourceFiles
    StantonNo.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_StantonNo_H
#define functionObjects_StantonNo_H

#include "fvMeshFunctionObject.H"
#include "writeFile.H"
#include "volFieldsFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                            Class StantonNo Declaration
\*---------------------------------------------------------------------------*/

class StantonNo
:
    public fvMeshFunctionObject,
    public writeFile
{
    // Private data

        //- Name of wallHeatFlux field, default is "wallHeatFlux"
        word wallHeatFlux_;
        
        //- Name of inflow patch
        word inflowPatchName_;
        
        //- header of wallHeatFlux field
        IOobject wallHeatFluxHeader_;
        
    
    // Private Member Functions

        //- File header information
        virtual void writeFileHeader(Ostream& os) const;

        //- Disallow default bitwise copy construct
        StantonNo(const StantonNo&);

        //- Disallow default bitwise assignment
        void operator=(const StantonNo&);


public:

    //- Runtime type information
    TypeName("StantonNo");


    // Constructors

        //- Construct from Time and dictionary
        StantonNo
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );


    //- Destructor
    virtual ~StantonNo();


    // Member Functions

        //- Read the StantonNo data
        virtual bool read(const dictionary&);

        //- Calculate the StantonNo field
        virtual bool execute();

        //- Write the StantonNo field
        virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
