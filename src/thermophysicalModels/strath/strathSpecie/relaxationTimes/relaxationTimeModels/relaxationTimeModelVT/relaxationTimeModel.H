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

Class
    Foam::relaxationTimeModel

Description
    Abstract base class for V-T relaxation models.

SourceFiles
    relaxationTimeModelI.H
    relaxationTimeModel.C

\*---------------------------------------------------------------------------*/

#ifndef relaxationTimeModel_H
#define relaxationTimeModel_H

#include "fvc.H"

#include "rho2ChemistryModel.H"
#include "turbulentFluidThermoModel.H"
#include "VTModel.H"
#include "multi2ComponentMixture.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class relaxationTimeModel Declaration
\*---------------------------------------------------------------------------*/

class relaxationTimeModel
:
    public IOdictionary
{

protected:

    // Protected data

        //- Mesh
        const fvMesh& mesh_;

        //- Reference to the thermodynamic package
        rho2ReactionThermo& thermo_;

        //- Reference to the turbulence model
        const compressibleTurbulenceModel& turbulence_;

        //- Interspecies V-T models
        autoPtr<VTModel> tauVTijModel_;

        //- V-T relaxation time [s]
        PtrList<volScalarField> tauVT_;
        //PtrList<PtrList<volScalarField> > tauVTmode_; // TODO ABORTIVE WORK

        //- V-T source term in the vib. energy equations [kg/m^3/s]
        PtrList<volScalarField> QVT_;
        //PtrList<PtrList<volScalarField> > QVTmode_; // TODO ABORTIVE WORK


private:

    // Private member functions

        //- Disallow copy construct
        relaxationTimeModel(const relaxationTimeModel&);

        //- Disallow default bitwise assignment
        void operator=(const relaxationTimeModel&);


public:

    //- Runtime type information
    TypeName("relaxationTimeModel");


    // Declare run-time constructor selection table

#ifndef SWIG
    declareRunTimeSelectionTable
    (
        autoPtr,
        relaxationTimeModel,
        fvMesh,
        (
            rho2ReactionThermo& thermo,
            const compressibleTurbulenceModel& turbulence
        ),
        (thermo, turbulence)
    );
#endif


    // Constructors

        //- Construct from components
        relaxationTimeModel
        (
            rho2ReactionThermo& thermo,
            const compressibleTurbulenceModel& turbulence
        );


    // Selectors

        //- Return a reference to the selected multispecies transport model

        static autoPtr<relaxationTimeModel> New
        (
            rho2ReactionThermo& thermo,
            const compressibleTurbulenceModel& turbulence
        );


    //- Destructor

        virtual ~relaxationTimeModel()
        {}


    // Member Functions

        // Access

            //- Return the const interspecies V-T relaxation time field
            inline const volScalarField& tauVTij(label i, label j) const;

            //- Return the const source term in the corresponding vib. energy
            //  equation due to V-T energy transfer
            inline const volScalarField& QVT(const label speciei) const;

            //- Return the const source term in the corresponding vib. energy
            //  equation due to V-T energy transfer
            inline const volScalarField& QVT(const word& specieName) const;

            //- Return the non-const source term in the corresponding vib.
            //  energy equation due to V-T energy transfer
            inline volScalarField& QVT(const label speciei);

            //- Return the non-const source term in the corresponding vib.
            //  energy equation due to V-T energy transfer
            inline volScalarField& QVT(const word& specieName);

            //- Return the const Source term in vib. energy equations
            //  due to V-T energy transfer, vibrational mode vibMode
            //inline const volScalarField& QVT_mode(label speciei, label vibMode) const;  // TODO ABORTIVE WORK

            //- Return the table of species
            inline const speciesTable& species() const;

            //- Return the table of species having a vib. eq. being solved
            //  or a hack of this table for the single-vibrational pool model
            // that will contain all (ionised/nonionised) molecules
            inline const speciesTable& solvedVibEqSpecies() const;

            //- Molar weight for a specie given by index
            virtual scalar W(label i) const = 0;


        //- Source term in vibro-electronic energy equation due to V-T energy
        //  transfer
        tmp<volScalarField> QVT();

        //- Correct the V-T relaxation times and V-T energy transfer
        virtual void correct() = 0;

        //- Read thermo2T dictionary
        virtual bool read() = 0;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "relaxationTimeModelI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
