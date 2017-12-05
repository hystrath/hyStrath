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
    polyStandardFields

Description
     These are 3D fields for measurements in MD simulations. 

    Cumulative averaging 
    
SourceFiles
    polyStandardFieldsI.H
    polyStandardFields.C

\*---------------------------------------------------------------------------*/

#ifndef polyStandardFields_H
#define polyStandardFields_H

#include "vector.H"
#include "fvMesh.H"
#include "polyMolecule.H"
#include "volFields.H"
#include "scalarIOField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class polyStandardFields Declaration
\*---------------------------------------------------------------------------*/
// class polyMoleculeCloud;


class polyStandardFields
{
    // Private data

        const fvMesh& mesh_;

        polyMoleculeCloud& cloud_;

        // field name
        word fieldName_;

        // ****
        List<label> molIds_;
        
        // geometric fields
        
        //- number density field
        volScalarField rhoN_;

        //- Mass density field
        volScalarField rhoM_;

        //- velocity field
        volVectorField USAM_;
        

        // collectors 
        scalarField mass_;
        vectorField momentum_;
        scalar nTimeSteps_;
        
//         volVectorField UCAM_;        

//         //- temperature
//         volScalarField T_;
// 
//         //- scalar pressure
//         volScalarField p_;
        


public:

//     friend class dsmcParcel;

    // Constructors

        //- Construct from mesh
        polyStandardFields
        (
            const fvMesh& mesh,
            polyMoleculeCloud& cloud,
            const word& fieldName,         
            bool write
        );

        //- Construct from mesh and cloud
        polyStandardFields
        (
            const fvMesh& mesh,
            polyMoleculeCloud& cloud,
            const word& fieldName
        );


    // Destructor

        ~polyStandardFields();

    // Member Functions

//         void calculateFields();

        void calculateFields();

        void resetFields();

        // Access

    // References to the boundary fields for surface data collection
/*
    //- Return non-const heat flux boundary field reference
    inline volScalarField::GeometricBoundaryField& qBF();

    //- Return non-const force density at boundary field reference
    inline volVectorField::GeometricBoundaryField& fDBF();

    //- Return non-const number density boundary field reference
    inline volScalarField::GeometricBoundaryField& rhoNBF();

    //- Return non-const mass density boundary field reference
    inline volScalarField::GeometricBoundaryField& rhoMBF();

    //- Return non-const linear kinetic energy density boundary
    //  field reference
    inline volScalarField::GeometricBoundaryField& linearKEBF();

    //- Return non-const rotational energy density boundary field
    // reference
    inline volScalarField::GeometricBoundaryField& rotationalEBF();

    //- Return non-const rotational degree of freedom density boundary
    //  field reference
    inline volScalarField::GeometricBoundaryField& rotationalDofBF();

    //- Return non-const vibrational energy density boundary field
    // reference
    inline volScalarField::GeometricBoundaryField& vibrationalEBF();

    //- Return non-const vibrational degree of freedom density boundary
    //  field reference
    inline volScalarField::GeometricBoundaryField& vibrationalDofBF();

    //- Return non-const momentum density boundary field reference
    inline volVectorField::GeometricBoundaryField& momentumBF();

    //- Return non-const boundaryT boundary field reference
    // 	inline volScalarField::GeometricBoundaryField& boundaryRadiationTBF();*/



    // References to the macroscopic fields


    //- Return heat flux at surface field
    inline const volScalarField& q() const;

    //- Return force density at surface field
    inline const volVectorField& fD() const;

    //- Return the real particle number density field
    inline const volScalarField& rhoN() const;

    //- Return the particle mass density field
    inline const volScalarField& rhoM() const;

    //- Return the field of number of DSMC particles
    inline const volScalarField& dsmcRhoN() const;

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

#include "polyStandardFieldsI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
