/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013-2016 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2015 OpenCFD Ltd.
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
    Foam::functionObjects::PecletNo

Group
    grpFieldFunctionObjects

Description
    Calculates and outputs the Peclet number as a surfaceScalarField.

Usage
    Example of function object specification to calculate the Peclet number:
    \verbatim
    PecletNo1
    {
        type        PecletNo;
        libs        ("libfieldFunctionObjects.so");
        ...
    }
    \endverbatim

    Where the entries comprise:
    \table
        Property     | Description             | Required    | Default value
        type         | type name: Peclet       | yes         |
        phi          | Name of flux field      | no          | phi
        rho          | Name of density field   | no          | rho
        result       | Name of Peclet field    | no          | <function name>
        log          | Log to standard output  | no          | yes
    \endtable

SourceFiles
    PecletNo.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_PecletNo_H
#define functionObjects_PecletNo_H

#include "fieldExpression.H"
#include "surfaceFieldsFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                          Class PecletNo Declaration
\*---------------------------------------------------------------------------*/

class PecletNo
:
    public fieldExpression
{
    // Private data

        //- Name of density field, default is "rho"
        word rhoName_;


    // Private Member Functions

        //- Optionally scale the flux for compressible cases
        tmp<surfaceScalarField> rhoScale(const surfaceScalarField& phi) const;

        //- Calculate the Peclet number field and return true if successful
        virtual bool calc();


public:

    //- Runtime type information
    TypeName("PecletNo");


    // Constructors

        //- Construct from Time and dictionary
        PecletNo
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );


    //- Destructor
    virtual ~PecletNo();


    // Member Functions

        //- Read the PecletNo data
        virtual bool read(const dictionary&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
