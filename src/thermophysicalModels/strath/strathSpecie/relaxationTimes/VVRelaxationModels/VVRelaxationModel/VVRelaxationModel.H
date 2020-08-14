/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016-2020 hyStrath
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of hyStrath, a derivative work of OpenFOAM.

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

Namespace
    Foam::VVRelaxationModels

Description
    Namespace for V-V relaxation models.

Class
    Foam::VVRelaxationModel

Description
    An abstract base class for V-V relaxation models.

SourceFiles
    VVRelaxationModel.C

\*---------------------------------------------------------------------------*/

#ifndef VVRelaxationModel_H
#define VVRelaxationModel_H

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
                Class VVRelaxationModel Declaration
\*---------------------------------------------------------------------------*/

class VVRelaxationModel
{

protected:

    // Protected data

        //- Two-temperature model dictionary
        const dictionary& dict1_;

        //- Thermophysical properties dictionary
        const dictionary& dict2_;

        //- Species name
        const word& name1_, name2_;

        //- Species label associated with thee names
        const label& lname1_, lname2_;

        // Pressure [Pa]
        const volScalarField& p_;

        // Translational temperature [K]
        const volScalarField& T_;

        // Vibrational temperatures [K]
        const PtrList<volScalarField>& Tv_;

        // Number densities [m^-3]
        const PtrList<volScalarField>& nD_;

        //- Overwrite the default coefficients for V-V energy transfers
        bool VVOverwriteDefault_;

        //- If overwrite is on, are the coefficients species dependent
        bool VVSpeciesDependent_;

        //- If overwrite is on, does the colliding pair matters
        bool VVCollidingPartner_;


    // Private Member Functions

        //- Disallow copy construct
        VVRelaxationModel(const VVRelaxationModel&);

        //- Disallow default bitwise assignment
        void operator=(const VVRelaxationModel&);


public:

    //- Runtime type information
    TypeName("VVRelaxationModel");


    // Declare run-time constructor selection table

      declareRunTimeSelectionTable
      (
          autoPtr,
          VVRelaxationModel,
          dictionary,
          (
              const word& name1,
              const word& name2,
              const label& lname1,
              const label& lname2,
              const dictionary& dict1,
              const dictionary& dict2,
              const volScalarField& p,
              const volScalarField& Tt,
              const PtrList<volScalarField>& Tv,
              const PtrList<volScalarField>& nD
          ),
          (name1, name2, lname1, lname2, dict1, dict2, p, Tt, Tv, nD)
      );


    // Selectors

        //- Return a reference to the selected V-V relaxation model
        static autoPtr<VVRelaxationModel> New
        (
            const word& name1,
            const word& name2,
            const label& lname1,
            const label& lname2,
            const dictionary& dict1,
            const dictionary& dict2,
            const volScalarField& p,
            const volScalarField& Tt,
            const PtrList<volScalarField>& Tv,
            const PtrList<volScalarField>& nD
        );


    // Constructors

        //- Construct from components
        VVRelaxationModel
        (
            const word& name1,
            const word& name2,
            const label& lname1,
            const label& lname2,
            const dictionary& dict1,
            const dictionary& dict2,
            const volScalarField& p,
            const volScalarField& Tt,
            const PtrList<volScalarField>& Tv,
            const PtrList<volScalarField>& nD
        );


    // Destructor

        virtual ~VVRelaxationModel()
        {}


    // Member Functions

        //- V-T relaxation time [s]
         virtual tmp<volScalarField> tauVV() const = 0;

        //- V-T relaxation time for patch [s]
        virtual tmp<scalarField> tauVV
        (
            const label patchi,
            const scalarField& p,
            const scalarField& Tt,
            const PtrList<scalarField>& Tv,
            const PtrList<scalarField>& nD
        ) const = 0;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
