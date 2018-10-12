/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2005 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Class
    pdStandardFields

Description

SourceFiles
    pdStandardFieldsI.H
    pdStandardFields.C

\*---------------------------------------------------------------------------*/

#ifndef pdStandardFields_H
#define pdStandardFields_H

#include "vector.H"
#include "fvMesh.H"
#include "pdParcel.H"
#include "volFields.H"
#include "scalarIOField.H"

#include "autoPtr.H"
#include "polyMeshTetDecomposition.H"
#include "AveragingMethod.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class pdStandardFields Declaration
\*---------------------------------------------------------------------------*/
// class pdCloud;


class pdStandardFields
{
    // Private data

        const fvMesh& mesh_;

        pdCloud& cloud_;

        //- setup averaging method
        autoPtr<AveragingMethod<scalar> > AveragePtr_;

        //- Heat flux at surface field
        volScalarField q_;

        //- Force density at surface field
        volVectorField fD_;

        //- number density field
        volScalarField rhoN_;

        //- Mass density field
        volScalarField rhoM_;

        //- pd particle density field
        volScalarField pdRhoN_;

        //- linear kinetic energy density field
        volScalarField linearKE_;

        //- Rotational energy density field
        volScalarField rotationalE_;

        // Rotational degree of freedom density field
        volScalarField rotationalDof_;

	//- Vibrational energy density field
        volScalarField vibrationalE_;

        // Vibrational degree of freedom density field
        volScalarField vibrationalDof_;

        //- Momentum density field
        volVectorField momentum_;

        //- Velocity field
        //volVectorField stdU_;

        //- translational temperature field
        //volScalarField transT_;

	//- boundaryT field
//         volScalarField boundaryRadiationT_;



public:

//     friend class pdParcel;

    // Constructors

        //- Construct from mesh
        pdStandardFields
        (
            const fvMesh& mesh,
            pdCloud& cloud,
            bool write
        );

        //- Construct from mesh and cloud
        pdStandardFields
        (
            const fvMesh& mesh,
            pdCloud& cloud
        );


    // Destructor

        ~pdStandardFields();

    // Member Functions

//         void calculateFields();

        void calculateFields();

        void resetFields();

        // Access

    // References to the boundary fields for surface data collection

    //- Return non-const heat flux boundary field reference
    inline volScalarField::Boundary& qBF();

    //- Return non-const force density at boundary field reference
    inline volVectorField::Boundary& fDBF();

    //- Return non-const number density boundary field reference
    inline volScalarField::Boundary& rhoNBF();

    //- Return non-const mass density boundary field reference
    inline volScalarField::Boundary& rhoMBF();

    //- Return non-const linear kinetic energy density boundary
    //  field reference
    inline volScalarField::Boundary& linearKEBF();

    //- Return non-const rotational energy density boundary field
    // reference
    inline volScalarField::Boundary& rotationalEBF();

    //- Return non-const rotational degree of freedom density boundary
    //  field reference
    inline volScalarField::Boundary& rotationalDofBF();

    //- Return non-const vibrational energy density boundary field
    // reference
    inline volScalarField::Boundary& vibrationalEBF();

    //- Return non-const vibrational degree of freedom density boundary
    //  field reference
    inline volScalarField::Boundary& vibrationalDofBF();

    //- Return non-const momentum density boundary field reference
    inline volVectorField::Boundary& momentumBF();

    //- Return non-const velocity boundary field reference
    //inline volVectorField::Boundary& stdUBF();

    //- Return non-const translational temperature boundary
    //  field reference
    //inline volScalarField::Boundary& transTBF();


    //- Return non-const boundaryT boundary field reference
    // 	inline volScalarField::Boundary& boundaryRadiationTBF();



    // References to the macroscopd fields


    //- Return heat flux at surface field
    inline const volScalarField& q() const;

    //- Return force density at surface field
    inline const volVectorField& fD() const;

    //- Return the real particle number density field
    inline const volScalarField& rhoN() const;

    //- Return the particle mass density field
    inline const volScalarField& rhoM() const;

    //- Return the field of number of PD particles
    inline const volScalarField& pdRhoN() const;

    //- Return the total linear kinetic energy (translational and
    // thermal density field
    inline const volScalarField& linearKE() const;

    //- Return the rotational energy density field
    inline const volScalarField& rotationalE() const;

    //- Return the average rotational degrees of freedom  field
    inline const volScalarField& rotationalDof() const;

    //- Return the vibrational energy density field
    inline const volScalarField& vibrationalE() const;

    //- Return the average vibrational degrees of freedom  field
    inline const volScalarField& vibrationalDof() const;

    //- Return the momentum density field
    inline const volVectorField& momentum() const;

    //- Return the velocity field
    //inline const volVectorField& stdU() const;

    //- Return the translational temperature  field
    //inline const volScalarField& transT() const;


	//- Return the boundary density field
// 	inline const volScalarField& boundaryRadiationT() const;


        // Check

        // Edit

        // Write

    // Member Operators

    // Friend Functions

    // Friend Operators


    // IOstream Operators

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "pdStandardFieldsI.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
