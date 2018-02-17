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
    Foam::mixed2EnergyFvPatchScalarField

Group
    grpThermoBoundaryConditions

Description
    This boundary condition provides a mixed condition for the internal energy

SeeAlso
    Foam::mixedFvPatchField

SourceFiles
    mixed2EnergyFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef mixed2EnergyFvPatchScalarField_H
#define mixed2EnergyFvPatchScalarField_H

#include "mixedFvPatchFields.H"

#include "rho2ReactionThermo.H" // NEW VINCENT 18/04/2016

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
               Class mixed2EnergyFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class mixed2EnergyFvPatchScalarField
:
    public mixedFvPatchScalarField
{

private:

    // Private data

        // Thermo package that contains species patch data
        const rho2ReactionThermo& thermo_;
        

public:

    //- Runtime type information
    TypeName("mixed2Energy");


    // Constructors

        //- Construct from patch and internal field
        mixed2EnergyFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        mixed2EnergyFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given mixed2EnergyFvPatchScalarField
        // onto a new patch
        mixed2EnergyFvPatchScalarField
        (
            const mixed2EnergyFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        mixed2EnergyFvPatchScalarField
        (
            const mixed2EnergyFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new mixed2EnergyFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        mixed2EnergyFvPatchScalarField
        (
            const mixed2EnergyFvPatchScalarField&,
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
                new mixed2EnergyFvPatchScalarField(*this, iF)
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
