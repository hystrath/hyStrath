/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2010 OpenCFD Ltd.
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
    Foam::nonLinMaxwellSlipUFvPatchVectorField

Description
    nonLinMaxwell slip boundary condition including thermal creep and surface
    curvature terms that can be optionally switched off.

SourceFiles
    fixedRhoFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef nonLinMaxwellSlipUFvPatchVectorField_H
#define nonLinMaxwellSlipUFvPatchVectorField_H

#include "mixedFixedValueSlipFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class nonLinMaxwellSlipUFvPatch Declaration
\*---------------------------------------------------------------------------*/

class nonLinMaxwellSlipUFvPatchVectorField
:
    public mixedFixedValueSlipFvPatchVectorField
{

    // Private data

        // Wall velocity
        vectorField Uwall_;

//        scalar refSlipLen_;
        scalar rho2_;
        scalar rho1_;
        scalar rho0_;
//
//        scalar gammaC5_;
//        scalar gammaC4_;
//        scalar gammaC3_;
        scalar gammaC2_;
        scalar gammaC1_;
//        scalar gammaC0_;

//        scalar p00_;
//        scalar p10_;
//        scalar p01_;
//
//        scalar p20_;
//        scalar p11_;
//        scalar p30_;
//        scalar p21_;

        // scalar frictionCoeff_;


public:

    //- Runtime type information
    TypeName("nonLinMaxwellSlipU");


    // Constructors

        //- Construct from patch and internal field
        nonLinMaxwellSlipUFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        nonLinMaxwellSlipUFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  nonLinMaxwellSlipUFvPatchVectorField onto a new patch
        nonLinMaxwellSlipUFvPatchVectorField
        (
            const nonLinMaxwellSlipUFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new nonLinMaxwellSlipUFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        nonLinMaxwellSlipUFvPatchVectorField
        (
            const nonLinMaxwellSlipUFvPatchVectorField&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchVectorField> clone
        (
            const DimensionedField<vector, volMesh>& iF
        ) const
        {
            return tmp<fvPatchVectorField>
            (
                new nonLinMaxwellSlipUFvPatchVectorField(*this, iF)
            );
        }

    // Member functions

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
