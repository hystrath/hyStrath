/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013-2016 OpenFOAM Foundation
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
    Foam::functionObjects::blendingFactor

Group
    grpFieldFunctionObjects

Description
    Calculates and outputs the blendingFactor as used by the bended convection
    schemes.  The output is a volume field (cells) whose value is calculated via
    the maximum blending factor for any cell face.

    The weight of a blended scheme is given by a function of the blending
    factor, f:

    \[f
        weight = f scheme1 + (1 - f) scheme2
    \f]

    The factor is a face-based quantity, which is converted to a cell-based
    quantity by assigning the minimum blending factor for any cell face.

    An indicator (volume) field, named <functionObjectName>:<fieldName>, is
    generated that is set to (1 - f), i.e. values of:
    - 0 represent scheme1 as active, and
    - 1 represent scheme2 as active.
    - intermediate values show the contribution to scheme2

    Additional reporting is written to the standard output, providing
    statistics as to the number of cells used by each scheme.

Usage
    Example of function object specification to calculate the blending factor:
    \verbatim
    blendingFactor1
    {
        type        blendingFactor;
        libs        ("libfieldFunctionObjects.so");

        ...

        // Name of field
        field       U;
    }
    \endverbatim

    Where the entries comprise:
    \table
        Property     | Description             | Required    | Default value
        type         | Type name: blendingFactor | yes       |
        phi          | Name of flux field      | no          | phi
        field        | Name of field to evaluate | yes       |
        tolerance    | Tolerance for number of blended cells | no | 0.001
        log          | Log to standard output  | no          | yes
    \endtable

See also
    Foam::functionObjects::fieldExpression
    Foam::functionObjects::fvMeshFunctionObject
    Foam::functionObjects::writeFile

SourceFiles
    blendingFactor.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_blendingFactor_H
#define functionObjects_blendingFactor_H

#include "fieldExpression.H"
#include "writeFile.H"
#include "convectionScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                       Class blendingFactor Declaration
\*---------------------------------------------------------------------------*/

class blendingFactor
:
    public fieldExpression,
    public writeFile
{
    // Private member data

        //- Name of flux field, default is "phi"
        word phiName_;

        //- Tolerance used when calculating the number of blended cells
        scalar tolerance_;


    // Private Member Functions

        //- Calculate the blending factor field
        template<class Type>
        void calcBlendingFactor
        (
            const GeometricField<Type, fvPatchField, volMesh>& field,
            const typename fv::convectionScheme<Type>& cs
        );

        //- Calculate the blending factor field
        template<class Type>
        bool calcScheme();

        //- Calculate the blending factor field and return true if successful
        virtual bool calc();


protected:

    // Protected Member Functions

        //- Write the file header
        virtual void writeFileHeader(Ostream& os) const;


public:

    //- Runtime type information
    TypeName("blendingFactor");


    // Constructors

        //- Construct from Time and dictionary
        blendingFactor
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );


    //- Destructor
    virtual ~blendingFactor();


    // Member Functions

        //- Read the blendingFactor data
        virtual bool read(const dictionary&);

        //- Write the blendingFactor
        virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "blendingFactorTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
