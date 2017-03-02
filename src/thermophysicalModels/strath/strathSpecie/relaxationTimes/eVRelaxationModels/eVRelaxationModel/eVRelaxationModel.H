/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
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

Namespace
    Foam::eVRelaxationModels

Description
    Namespace for e-V relaxation models.

Class
    Foam::eVRelaxationModel

Description
    An abstract base class for e-V relaxation models.

SourceFiles
    eVRelaxationModel.C

\*---------------------------------------------------------------------------*/

#ifndef eVRelaxationModel_H
#define eVRelaxationModel_H

#include "advancedSpecie.H"
#include "IOdictionary.H"
#include "typeInfo.H"
#include "runTimeSelectionTables.H"
#include "scalarField.H"
#include "volFieldsFwd.H"
#include "surfaceFieldsFwd.H"
#include "dimensionedScalar.H"
#include "tmp.H"
#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class eVRelaxationModel Declaration
\*---------------------------------------------------------------------------*/

class eVRelaxationModel
{

protected:

    // Protected data
    
        //- Two-temperature model dictionary
        const dictionary& dict1_;
        
        //- Thermophysical properties dictionary
        const dictionary& dict2_;

        //- Species name
        const word& name1_;
        
        //- Species label associated with the names
        const label& lname1_;

        // Electron pressure [Pa]
        const volScalarField& p_;

        // Vibrational temperature [K]
        const volScalarField& Tv_;
        
        //- Overwrite the default coefficients for e-V energy transfers
        bool eVOverwriteDefault_;
        
        //- If overwrite is on, are the coefficients species dependent
        bool eVSpeciesDependent_;
        

    // Private Member Functions

        //- Disallow copy construct
        eVRelaxationModel(const eVRelaxationModel&);

        //- Disallow default bitwise assignment
        void operator=(const eVRelaxationModel&);


public:

    //- Runtime type information
    TypeName("eVRelaxationModel");


    // Declare run-time constructor selection table
    
      declareRunTimeSelectionTable
      (
          autoPtr,
          eVRelaxationModel,
          dictionary,
          (
              const word& name1,
              const label& lname1,
              const dictionary& dict1,
              const dictionary& dict2,
              const volScalarField& p,
              const volScalarField& Tv
          ),
          (name1, lname1, dict1, dict2, p, Tv)
      );


    // Selectors

        //- Return a reference to the selected e-V relaxation model
        static autoPtr<eVRelaxationModel> New
        (
            const word& name1,
            const label& lname1,
            const dictionary& dict1,
            const dictionary& dict2,
            const volScalarField& p,
            const volScalarField& Tv
        );


    // Constructors

        //- Construct from components
        eVRelaxationModel
        (
            const word& name1,
            const label& lname1,
            const dictionary& dict1,
            const dictionary& dict2,
            const volScalarField& p,
            const volScalarField& Tv
        );


    // Destructor

        virtual ~eVRelaxationModel()
        {}


    // Member Functions
    
        //- e-V relaxation time [s]
         virtual tmp<volScalarField> taueV() const = 0;
         
        //- e-V relaxation time for patch [s]
        virtual tmp<scalarField> taueV
        (
            const label patchi,
            const scalarField& p,
            const scalarField& Tv
        ) const = 0;
       
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
