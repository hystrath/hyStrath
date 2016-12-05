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
    Foam::mixed2VELEnergyFvPatchScalarField

Group
    grpThermoBoundaryConditions

Description
    This boundary condition provides a mixed condition for the 
    vibro-electronic mode of the internal energy

SeeAlso
    Foam::mixedFvPatchField

SourceFiles
    mixed2VELEnergyFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef mixed2VELEnergyFvPatchScalarField_H
#define mixed2VELEnergyFvPatchScalarField_H

#include "mixedFvPatchFields.H"

#include "rho2ReactionThermo.H" // NEW VINCENT 18/04/2016

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
               Class mixed2VELEnergyFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class mixed2VELEnergyFvPatchScalarField
:
    public mixedFvPatchScalarField
{

private:

    // Private data

        //- Species index
        word specieName_;
        
        // Thermo package that contains species patch data
        const rho2ReactionThermo& thermo_;
        

public:

    //- Runtime type information
    TypeName("mixed2VELEnergy");


    // Constructors

        //- Construct from patch and internal field
        mixed2VELEnergyFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        mixed2VELEnergyFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given mixed2VELEnergyFvPatchScalarField
        // onto a new patch
        mixed2VELEnergyFvPatchScalarField
        (
            const mixed2VELEnergyFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        mixed2VELEnergyFvPatchScalarField
        (
            const mixed2VELEnergyFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new mixed2VELEnergyFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        mixed2VELEnergyFvPatchScalarField
        (
            const mixed2VELEnergyFvPatchScalarField&,
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
                new mixed2VELEnergyFvPatchScalarField(*this, iF)
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
