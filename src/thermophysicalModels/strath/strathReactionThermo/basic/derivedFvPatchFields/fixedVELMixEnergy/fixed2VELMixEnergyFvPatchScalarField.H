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
    Foam::fixed2VELMixEnergyFvPatchScalarField

Group
    grpThermoBoundaryConditions

Description
    This boundary condition provides a fixed condition for internal energy

    \heading Patch usage

    Example of the boundary condition specification:
    \verbatim
    myPatch
    {
        type            fixed2VELMixEnergy;
        value           uniform 100;
    }
    \endverbatim

SeeAlso
    Foam::fixedValueFvPatchField

SourceFiles
    fixed2VELMixEnergyFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef fixed2VELMixEnergyFvPatchScalarField_H
#define fixed2VELMixEnergyFvPatchScalarField_H

#include "fixedValueFvPatchFields.H"

#include "rho2ReactionThermo.H" // NEW VINCENT 15/02/2017

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
               Class fixed2VELMixEnergyFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class fixed2VELMixEnergyFvPatchScalarField
:
    public fixedValueFvPatchScalarField
{

private:

    // Private data

        // Thermo package that contains species patch data
        const rho2ReactionThermo& thermo_; // NEW VINCENT 15/02/2017

public:

    //- Runtime type information
    TypeName("fixed2VELMixEnergy");


    // Constructors

        //- Construct from patch and internal field
        fixed2VELMixEnergyFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        fixed2VELMixEnergyFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given fixed2VELMixEnergyFvPatchScalarField
        // onto a new patch
        fixed2VELMixEnergyFvPatchScalarField
        (
            const fixed2VELMixEnergyFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        fixed2VELMixEnergyFvPatchScalarField
        (
            const fixed2VELMixEnergyFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new fixed2VELMixEnergyFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        fixed2VELMixEnergyFvPatchScalarField
        (
            const fixed2VELMixEnergyFvPatchScalarField&,
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
                new fixed2VELMixEnergyFvPatchScalarField(*this, iF)
            );
        }


    // Member functions

        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


#endif

// ************************************************************************* //
