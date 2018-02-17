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
    Foam::outletMappedUniformInletHeatAdditionFvPatchField

Group
    grpInletBoundaryConditions

Description

    This temperature boundary condition averages the temperature over the
    "outlet" patch specified by name "outletPatchName" and applies an extra
    heat source. This is set as a uniform temperature value on this patch.
    Additionally TMin/TMax limits can be applied

    \heading Patch usage

    \table
        Property    | Description             | Required    | Default value
        outletPatch | name of outlet patch    | yes         |
        Q           | Heat addition           | yes
        phi         | flux field name         | no          | phi
        TMin        | min temperature limit   | no          | 0.0
        TMax        | max temperature limit   | no          | 5000.0
    \endtable

    Example of the boundary condition specification:
    \verbatim
    myPatch
    {
        type            outletMappedUniformInletHeatAddition;
        outletPatch     aPatch;
        Q               10;     // Heat addition in W
        phi             phi;
        value           uniform 0;
    }
    \endverbatim

SeeAlso
    Foam::fixedValueFvPatchField

SourceFiles
    outletMappedUniformInletHeatAdditionFvPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef outletMappedUniformInletHeatAdditionFvPatchField_H
#define outletMappedUniformInletHeatAdditionFvPatchField_H

#include "fixedValueFvPatchFields.H"
#include "fvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
       Class outletMappedUniformInletHeatAdditionFvPatchField Declaration
\*---------------------------------------------------------------------------*/

class outletMappedUniformInletHeatAdditionFvPatchField
:
    public fixedValueFvPatchScalarField
{
    // Private data

        //- Name of the outlet patch to be mapped
        word outletPatchName_;

        //- Name of the flux transporting the field
        word phiName_;

        //- Input energy
        scalar Q_;

        //- Minimum Temperature Limit
        scalar TMin_;

        //- Maxmum Temperature Limit
        scalar TMax_;


public:

    //- Runtime type information
    TypeName("outletMappedUniformInletHeatAddition");


    // Constructors

        //- Construct from patch and internal field
        outletMappedUniformInletHeatAdditionFvPatchField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        outletMappedUniformInletHeatAdditionFvPatchField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given outletMappedUniformInletHeatAdditionFvPatchField
        //  onto a new patch
        outletMappedUniformInletHeatAdditionFvPatchField
        (
            const outletMappedUniformInletHeatAdditionFvPatchField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        outletMappedUniformInletHeatAdditionFvPatchField
        (
            const outletMappedUniformInletHeatAdditionFvPatchField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchField<scalar>> clone() const
        {
            return tmp<fvPatchField<scalar>>
            (
                new outletMappedUniformInletHeatAdditionFvPatchField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        outletMappedUniformInletHeatAdditionFvPatchField
        (
            const outletMappedUniformInletHeatAdditionFvPatchField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchField<scalar>> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchField<scalar>>
            (
                new outletMappedUniformInletHeatAdditionFvPatchField(*this, iF)
            );
        }


    // Member functions

        // Access

            //- Name of the outlet patch to be mapped
            const word& outletPatchName() const
            {
                return outletPatchName_;
            }


        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();


        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
