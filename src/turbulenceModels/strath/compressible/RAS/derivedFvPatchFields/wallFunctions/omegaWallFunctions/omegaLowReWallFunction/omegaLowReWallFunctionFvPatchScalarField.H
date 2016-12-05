/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012 OpenFOAM Foundation
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
    Foam::compressible::omegaLowReWallFunctionFvPatchScalarField

Group
    grpCmpWallFunctions

Description
    This boundary condition provides a turbulence kinetic energy wall function
    condition for low- and high-Reynolds number turbulent flow cases.

    The model operates in two modes, based on the computed laminar-to-turbulent
    switch-over y+ value derived from kappa and E.

    \heading Patch usage

    \table
        Property     | Description             | Required    | Default value
        beta1        | model coefficient       | no          | 0.075
    \endtable

    Example of the boundary condition specification:
    \verbatim
    myPatch
    {
        type            omegaLowReWallFunction;
    }
    \endverbatim

SeeAlso
    Foam::fixedValueFvPatchField

SourceFiles
    omegaLowReWallFunctionFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef omegaLowReWallFunctionFvPatchScalarField_H
#define omegaLowReWallFunctionFvPatchScalarField_H

#include "fixedValueFvPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{

/*---------------------------------------------------------------------------*\
            Class omegaLowReWallFunctionFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class omegaLowReWallFunctionFvPatchScalarField
:
    public fixedValueFvPatchField<scalar>
{
protected:

    // Protected data

        //- beta1 coefficient
        scalar beta1_;
        
        //- Von Karman constant
        scalar kappa_;

        //- E coefficient
        scalar E_;

        //- Y+ at the edge of the laminar sublayer
        scalar yPlusLam_;


    // Protected Member Functions

        //- Check the type of the patch
        virtual void checkType();

        //- Calculate the Y+ at the edge of the laminar sublayer
        scalar yPlusLam(const scalar kappa, const scalar E);


public:

    //- Runtime type information
    TypeName("compressible::omegaLowReWallFunction");


    // Constructors

        //- Construct from patch and internal field
        omegaLowReWallFunctionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        omegaLowReWallFunctionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given omegaLowReWallFunctionFvPatchScalarField
        //  onto a new patch
        omegaLowReWallFunctionFvPatchScalarField
        (
            const omegaLowReWallFunctionFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        omegaLowReWallFunctionFvPatchScalarField
        (
            const omegaLowReWallFunctionFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new omegaLowReWallFunctionFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        omegaLowReWallFunctionFvPatchScalarField
        (
            const omegaLowReWallFunctionFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new omegaLowReWallFunctionFvPatchScalarField(*this, iF)
            );
        }


    // Member functions

        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();

            //- Evaluate the patchField
            virtual void evaluate(const Pstream::commsTypes);


        // I-O

            //- Write
            virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace compressible
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
