/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 3 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Class
    Foam::radiation::fixedIncidentRadiationFvPatchScalarField

Group
    grpThermoBoundaryConditions

Description
    Boundary condition for thermal coupling for solid regions.
    Used to emulate a fixed incident radiative heat flux on a wall.

    the gradient heat flux is calculated as :

        Qr = emissivity*(QrIncident - sigma_*T^4)

    where:

    emissivity is the emissivity of the solid.
    QrIncident is the specified fixed incident radiation.

    Example usage:

    wall
    {
        type            fixedIncidentRadiation;
        QrIncident      uniform 500;
        kappa           solidThermo;
        KappaName       none;
    }

    kappa:
    - 'lookup' :    lookup volScalarField (or volSymmTensorField) with name
    - 'solidThermo' : use solidThermo kappa()

    emissivity:
    - 'lookup' :     lookup volScalarField emissivity
    - 'localSolidRadiation': Look up for local solidRadiation


SourceFiles
    fixedIncidentRadiationFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef fixedIncidentRadiationFvPatchScalarField_H
#define fixedIncidentRadiationFvPatchScalarField_H

#include "fixedGradientFvPatchFields.H"
#include "temperatureCoupledBase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace radiation
{
/*---------------------------------------------------------------------------*\
        Class fixedIncidentRadiationFvPatchScalarField declaration
\*---------------------------------------------------------------------------*/

class fixedIncidentRadiationFvPatchScalarField
:
    public fixedGradientFvPatchScalarField,
    public temperatureCoupledBase
{
    // Private data

        //- Incident radiative heat flux
        scalarField QrIncident_;


public:

    //- Runtime type information
    TypeName("fixedIncidentRadiation");


    // Constructors

        //- Construct from patch and internal field
        fixedIncidentRadiationFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        fixedIncidentRadiationFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  turbulentTemperatureCoupledBaffleMixedFvPatchScalarField onto a
        //  new patch
        fixedIncidentRadiationFvPatchScalarField
        (
            const
            fixedIncidentRadiationFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        fixedIncidentRadiationFvPatchScalarField
        (
            const fixedIncidentRadiationFvPatchScalarField&
        );


        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new fixedIncidentRadiationFvPatchScalarField
                (
                    *this
                )
            );
        }

        //- Construct as copy setting internal field reference
        fixedIncidentRadiationFvPatchScalarField
        (
            const fixedIncidentRadiationFvPatchScalarField&,
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
                new fixedIncidentRadiationFvPatchScalarField
                (
                    *this,
                    iF
                )
            );
        }


    // Member functions


     // Mapping functions

        //- Map (and resize as needed) from self given a mapping object
        virtual void autoMap(const fvPatchFieldMapper&);

        //- Reverse map the given fvPatchField onto this fvPatchField
        virtual void rmap
        (
            const fvPatchScalarField&,
            const labelList&
        );

        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam
} // End namespace radiation

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
