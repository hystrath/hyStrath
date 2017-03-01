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
    Foam::mixed2VELMixEnergyFvPatchScalarField

Group
    grpThermoBoundaryConditions

Description
    This boundary condition provides a mixed condition for the 
    trnas-rotational mode of the internal energy

SeeAlso
    Foam::mixedFvPatchField

SourceFiles
    mixed2VELMixEnergyFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef mixed2VELMixEnergyFvPatchScalarField_H
#define mixed2VELMixEnergyFvPatchScalarField_H

#include "mixedFvPatchFields.H"

#include "rho2ReactionThermo.H" // NEW VINCENT 15/02/2017

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
               Class mixed2VELMixEnergyFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class mixed2VELMixEnergyFvPatchScalarField
:
    public mixedFvPatchScalarField
{

private:

    // Private data

        // Thermo package that contains species patch data
        const rho2ReactionThermo& thermo_; // NEW VINCENT 15/02/201
        

public:

    //- Runtime type information
    TypeName("mixed2VELMixEnergy");


    // Constructors

        //- Construct from patch and internal field
        mixed2VELMixEnergyFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        mixed2VELMixEnergyFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given mixed2VELMixEnergyFvPatchScalarField
        // onto a new patch
        mixed2VELMixEnergyFvPatchScalarField
        (
            const mixed2VELMixEnergyFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        mixed2VELMixEnergyFvPatchScalarField
        (
            const mixed2VELMixEnergyFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new mixed2VELMixEnergyFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        mixed2VELMixEnergyFvPatchScalarField
        (
            const mixed2VELMixEnergyFvPatchScalarField&,
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
                new mixed2VELMixEnergyFvPatchScalarField(*this, iF)
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
