/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
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
    Foam::fvConductorPatch

SourceFiles
    fvConductorPatch.C

SeeAlso
    Foam::fixedValueFvPatchField

\*---------------------------------------------------------------------------*/

#ifndef fvConductorPatch_H
#define fvConductorPatch_H

#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
             Class fvConductorPatch Declaration
\*---------------------------------------------------------------------------*/

class fvConductorPatch
:
    public fixedValueFvPatchScalarField
{
    // Private data

        //- Name of the potentialField field
        word phiEName_;

        word wallQName_;

        scalar r_;

        scalar epsilonM_;

        scalarField phiW_;


public:

    //- Runtime type information
    TypeName("fvConductorPatch");


    // Constructors

        //- Construct from patch and internal field
        fvConductorPatch
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        fvConductorPatch
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given fvConductorPatch
        //  onto a new patch
        fvConductorPatch
        (
            const fvConductorPatch&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        fvConductorPatch
        (
            const fvConductorPatch&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new fvConductorPatch(*this)
            );
        }

        //- Construct as copy setting internal field reference
        fvConductorPatch
        (
            const fvConductorPatch&,
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
                new fvConductorPatch(*this, iF)
            );
        }


    // Member functions

        // Access

            //- Return the total pressure
            const scalarField& phiW() const
            {
                return phiW_;
            }

            //- Return reference to the total pressure to allow adjustment
            scalarField& phiW()
            {
                return phiW_;
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
                const fvPatchScalarField&,
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
