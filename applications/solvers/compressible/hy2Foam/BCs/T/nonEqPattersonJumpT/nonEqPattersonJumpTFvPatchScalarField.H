/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016-2020 hyStrath
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of hyStrath, a derivative work of OpenFOAM.

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
    Foam::nonEqPattersonJumpTFvPatchScalarField

Description
    Patterson temperature jump boundary condition

SourceFiles
    nonEqPattersonJumpTFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef nonEqPattersonJumpTFvPatchScalarFields_H
#define nonEqPattersonJumpTFvPatchScalarFields_H

#include "mixedFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class smoluchowskiJumpTFvPatch Declaration
\*---------------------------------------------------------------------------*/

class nonEqPattersonJumpTFvPatchScalarField
:
    public mixedFvPatchScalarField
{

    // Private data

        //- Velocity field name, default = "U"
        word UName_;

        //- Density field name, default = "rho"
        word rhoName_;

        //- Trans-rotational temperature, default = "Tt"
        word TName_;

        //- Dynamic viscosity field name, default = "mu"
        word muName_;

        //- Thermal diffusivity field name, default = "alphatr"
        word alphaName_; // NEW VINCENT 03/03/2016

        //- Ratio of specific heat capacities field name, default = "gammatr"
        word gammatrName_; // NEW VINCENT 03/03/2016

        //- Mean free path field name, default = "mfp"
        word mfpName_; // NEW VINCENT 28/02/2016

        //- Accommodation coefficient
        scalar accommodationCoeff_;

        //- Wall surface temperature
        scalarField Twall_;

        //- Heat capacity ratio (default 1.4)
        //scalar gamma_;

public:

    //- Runtime type information
    TypeName("nonEqPattersonJumpT");


    // Constructors

        //- Construct from patch and internal field
        nonEqPattersonJumpTFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        nonEqPattersonJumpTFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given nonEqPattersonJumpTFvPatchScalarField
        //  onto a new patch
        nonEqPattersonJumpTFvPatchScalarField
        (
            const nonEqPattersonJumpTFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new nonEqPattersonJumpTFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        nonEqPattersonJumpTFvPatchScalarField
        (
            const nonEqPattersonJumpTFvPatchScalarField&,
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
                new nonEqPattersonJumpTFvPatchScalarField(*this, iF)
            );
        }

        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap
            (
                const fvPatchFieldMapper&
            );

            //- Reverse map the given fvPatchField onto this fvPatchField
            virtual void rmap
            (
                const fvPatchField<scalar>&,
                const labelList&
            );


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
