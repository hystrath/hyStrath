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
    Foam::VTRelaxationModels

Description
    Namespace for V-T relaxation models.

Class
    Foam::VTRelaxationModel

Description
    An abstract base class for inter-species V-T relaxation models.

SourceFiles
    VTRelaxationModel.C

\*---------------------------------------------------------------------------*/

#ifndef VTRelaxationModel_H
#define VTRelaxationModel_H

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
                Class VTRelaxationModel Declaration
\*---------------------------------------------------------------------------*/

class VTRelaxationModel
{

protected:

    // Protected data

        //- Two-temperature model dictionary
        const dictionary& dict1_;

        //- Thermophysical properties dictionary
        const dictionary& dict2_;

        //- Species name
        const word& name1_, name2_;

        //- Species label associated with the names
        const label& lname1_, lname2_;

        // Pressure [Pa]
        const volScalarField& p_;

        // Heavy-particle temperature [K]
        const volScalarField& T_;

        // Vibrational temperatures [K]
        const PtrList<volScalarField>& Tv_;

        // Number densities [m^-3]
        const PtrList<volScalarField>& nD_;

        //- Full (true) or partial/pre coefficients formulation
        bool VTFullCoeffsForm_;

        //- Overwrite the default coefficients for V-T energy transfers
        bool VTOverwriteDefault_;

        //- If overwrite is on, are the coefficients species dependent
        bool VTSpeciesDependent_;

        //- If overwrite is on, does the colliding pair matters
        bool VTCollidingPartner_;


    // Private Member Functions

        //- Disallow copy construct
        VTRelaxationModel(const VTRelaxationModel&);

        //- Disallow default bitwise assignment
        void operator=(const VTRelaxationModel&);


public:

    //- Runtime type information
    TypeName("VTRelaxationModel");


    // Declare run-time constructor selection table

      declareRunTimeSelectionTable
      (
          autoPtr,
          VTRelaxationModel,
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

        //- Return a reference to the selected V-T relaxation model
        static autoPtr<VTRelaxationModel> New
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
        VTRelaxationModel
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

        virtual ~VTRelaxationModel()
        {}


    // Member Functions

        //- V-T relaxation time [s]
         virtual tmp<volScalarField> tauVT() const = 0;

        //- V-T relaxation time for patch [s]
        virtual tmp<scalarField> tauVT
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
