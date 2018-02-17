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
    Foam::gradient2EnergyFvPatchScalarField

Group
    grpThermoBoundaryConditions

Description
    This boundary condition provides a gradient condition for internal energy,
    where the gradient is calculated using:

        \f[
            \nabla(e_p) = \nabla_\perp C_p(p, T) + \frac{e_p - e_c}{\Delta}
        \f]

    where
    \vartable
        e_p     | energy at patch faces [J]
        e_c     | energy at patch internal cells [J]
        p       | pressure [bar]
        T       | temperature [K]
        C_p     | specific heat [J/kg/K]
        \Delta  | distance between patch face and internal cell centres [m]
    \endvartable

    \heading Patch usage

    Example of the boundary condition specification:
    \verbatim
    myPatch
    {
        type            gradient2Energy;
        gradient        uniform 10;
    }
    \endverbatim

SeeAlso
    Foam::fixedGradientFvPatchField

SourceFiles
    gradient2EnergyFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef gradient2EnergyFvPatchScalarField_H
#define gradient2EnergyFvPatchScalarField_H

#include "fixedGradientFvPatchFields.H"

#include "rho2ReactionThermo.H" // NEW VINCENT 22/02/2016

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
             Class gradient2EnergyFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class gradient2EnergyFvPatchScalarField
:
    public fixedGradientFvPatchScalarField
{

private:

    // Private data

        // Thermo package that contains species patch data
        const rho2ReactionThermo& thermo_;


public:

    //- Runtime type information
    TypeName("gradient2Energy");


    // Constructors

        //- Construct from patch and internal field
        gradient2EnergyFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        gradient2EnergyFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given gradient2EnergyFvPatchScalarField
        // onto a new patch
        gradient2EnergyFvPatchScalarField
        (
            const gradient2EnergyFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        gradient2EnergyFvPatchScalarField
        (
            const gradient2EnergyFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new gradient2EnergyFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        gradient2EnergyFvPatchScalarField
        (
            const gradient2EnergyFvPatchScalarField&,
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
                new gradient2EnergyFvPatchScalarField(*this, iF)
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
