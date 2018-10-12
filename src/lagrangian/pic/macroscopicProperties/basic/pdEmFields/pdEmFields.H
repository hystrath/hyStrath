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
    pdEmFields

Description

SourceFiles
    pdEmFieldsI.H
    pdEmFields.C

\*---------------------------------------------------------------------------*/

#ifndef pdEmFields_H
#define pdEmFields_H

//#include "vector.H"
#include "fvMesh.H"
#include "pdParcel.H"

#include "AveragingMethod.H"
//#include "nBasic.H"
//#include "nDual.H"

#include "pointMesh.H"
#include "pointPatchField.H"
#include "autoPtr.H"

#include "volFields.H"
#include "scalarIOField.H"
#include "IOdictionary.H"

#include "polyMeshTetDecomposition.H"

#include "fvCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class pdEmFields Declaration
\*---------------------------------------------------------------------------*/
// class pdCloud;

//template<class Type>
//class nAveragingMethod;

class pdEmFields
{
    // Private data

        //- mesh referance
        const fvMesh& mesh_;

        //- cloud referance
        pdCloud& cloud_;

        //- electric field average pointer
        autoPtr<AveragingMethod<vector>> EAveragePtr_;

        //- magnetic field average pointer
        autoPtr<AveragingMethod<vector>> BAveragePtr_;

        //- current density field average pointer
        //autoPtr<AveragingMethod<vector> > JpAveragePtr_;

        //- current potential field average pointer
        autoPtr<AveragingMethod<scalar>> phiEAveragePtr_;

        //- charge density average pointer
        autoPtr<AveragingMethod<scalar>> rhoAveragePtr_;

public:

    //public data
        //- electric field V/m
        volVectorField E_;

        //- magnetic field N/(m * A)
        volVectorField B_;

        //- Electric potential field V
        volScalarField phiE_;

        //- particle current density field A/m^2
        volVectorField Jp_;

        //- total charge density field C/m^3
        volScalarField rhoQ_;

        //- total charge on wall patches C
        volScalarField wallQ_;

        //- electron charge fluid density field C/m^3
        volScalarField rhoEF_;

        //- debye Length
        //volScalarField lambdaD_;

        //- debug switch for charge
        bool debugCharge_;

        //-debug switch for fields
        bool debugField_;

//     friend class pdParcel;

    // Constructors

        //- Construct from mesh
        pdEmFields
        (
            const fvMesh& mesh,
            pdCloud& cloud,
            bool write
        );

        //- Construct from mesh and cloud
        pdEmFields
        (
            const fvMesh& mesh,
            pdCloud& cloud
        );


    // Destructor

        ~pdEmFields();

    // Member Functions

        //Calculate electromagnetic fields based on particle distribution
        void calculateFields();

        //Calculate forces on particles based on fields
        void calculateForces();

        //Sets up initialised particles to leapfrog intergration
        void setupLeapfrog();

        //Reset fields
        void resetFields();

        //void updateAverages();

        // Access

    // References to the boundary fields for surface data collection

    //- Return non-const charge density at boundary field reference
    //inline volScalarField::Boundary& lambdaDBF();

    //- Return non-const electric potential boundary field reference
    inline volScalarField::Boundary& rhoQBF();

    //- Return non-const electric potential boundary field reference
    inline volScalarField::Boundary& wallQBF();

    //- Return non-const electon fluid boundary field reference
    inline volScalarField::Boundary& rhoEFBF();

    //- Return non-const electric potential boundary field reference
    inline volScalarField::Boundary& phiEBF();

    //- Return non-const electric field at boundary field reference
    inline volVectorField::Boundary& EBF();

        //- Return non-const electric field at boundary field reference
    inline volVectorField::Boundary& BBF();

    //- Return non-const current density field at boundary field reference
    inline volVectorField::Boundary& JpBF();



    // References to the macroscopd fields

    //- Return current density at surface field
    //inline const volScalarField& lambdaD() const;

    //- Return total charge density field
    inline const volScalarField& rhoQ() const;

    //- Return total charge density field
    inline const volScalarField& wallQ() const;

    //- Return electron fluid charge density field
    inline const volScalarField& rhoEF() const;

    //- Return electric potential field
    inline const volScalarField& phiE() const;

    //- Return electric field
    inline const volVectorField& E() const;

    //- Return magnetic field
    inline const volVectorField& B() const;

    //- Return current density field
    inline const volVectorField& Jp() const;


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

#include "pdEmFieldsI.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
