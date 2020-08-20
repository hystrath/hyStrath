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
    Foam::multiSpeciesTransportModel

Description
    Abstract base class for multispecies mass transport model models.

SourceFiles
    multiSpeciesTransportModelI.H
    multiSpeciesTransportModel.C

\*---------------------------------------------------------------------------*/

#ifndef multiSpeciesTransportModel_H
#define multiSpeciesTransportModel_H

#include "fvc.H"

#include "rho2ChemistryModel.H"
#include "turbulentFluidThermoModel.H"
#include "diffusionModel.H"
#include "multi2ComponentMixture.H"

#include "basic2MultiComponentMixture.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class multiSpeciesTransportModel Declaration
\*---------------------------------------------------------------------------*/

class multiSpeciesTransportModel
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

        //- Species mass fluxes [kg/s]
        PtrList<surfaceScalarField> spMassFlux_;

        //- diffusion models
        autoPtr<diffusionModel> DijModel_;

        //- Multispecies diffusion coefficients [kg/m/s]
        PtrList<volScalarField> D_;

        //- Species diffusion fluxes [kg/m^2/s]
        PtrList<volVectorField> JnonCorrected_;

        //- Sum of the diffusion fluxes to enter into the species and reaction
        //  equation as an explicit divergence term [kg/m^2/s]
        volVectorField sumDiffusionFluxes_;

        //- Species diffusion fluxes due to the gradient of pressure [kg/m^2/s]
        PtrList<volVectorField> JGradp_;

        //- Species diffusion fluxes due to the gradient of temperature
        // [kg/m^2/s]
        PtrList<volVectorField> JGradT_;

        //- Use the non-corrected version of the mass diffusion term
        bool useNonCorrected_;

        //- Use the gradient of molar fraction in the species and transport
        //  equation and use the gradient of mass fraction if disabled
        //  Always true for the SCEBD model
        bool solvingForX_;

        //- Account of the effect of the gradient of pressure on mass diffusion
        bool addPressureGradientTerm_;

        //- Account of the effect of the gradient of temperature on mass
        //  diffusion (Soret effect)
        bool addTemperatureGradientTerm_;


    // Protected member functions

        //- Calculates and stores the non-corrected diffusion mass flux field
        //  for a specie given by index
        virtual void calculateJ(const label i);

        //- Summation of the diffusion fluxes to later ensure mass conservation
        virtual void calculateSumDiffusionFluxes();

        //- Return the corrected diffusion mass flux field for a specie given
        //  by index
        virtual volVectorField Jcorrected(const label i) const;

        //- Calculates the contribution of the pressure gradient for each
        //  species and stores it into JGradp_
        void pressureGradientContributionToSpeciesMassFlux();

        //- Calculates the contribution of the temperature gradient for each
        //  species and stores it into JGradT_
        void temperatureGradientContributionToSpeciesMassFlux();

private:

    // Private member functions

        //- Disallow copy construct
        multiSpeciesTransportModel(const multiSpeciesTransportModel&);

        //- Disallow default bitwise assignment
        void operator=(const multiSpeciesTransportModel&);


public:

    //- Runtime type information
    TypeName("multiSpeciesTransportModel");


    // Declare run-time constructor selection table

#ifndef SWIG
    declareRunTimeSelectionTable
    (
        autoPtr,
        multiSpeciesTransportModel,
        fvMesh,
        (
            rho2ReactionThermo& thermo,
            const compressible::turbulenceModel& turbulence
        ),
        (thermo, turbulence)
    );
#endif


    // Constructors

        //- Construct from components
        multiSpeciesTransportModel
        (
            rho2ReactionThermo& thermo,
            const compressible::turbulenceModel& turbulence
        );


    // Selectors

        //- Return a reference to the selected multispecies transport model

        static autoPtr<multiSpeciesTransportModel> New
        (
            rho2ReactionThermo& thermo,
            const compressible::turbulenceModel& turbulence
        );


    //- Destructor

        virtual ~multiSpeciesTransportModel()
        {}


    // Member Functions

        //- Term in the total energy equation due to diffusion fluxes
        volVectorField multiSpeciesHeatSource() const;

        //- Update the species diffusion coefficients
        virtual void correct() = 0;

        //- Get the mass flux for specie i and store it into spMassFlux
        void getSpeciesMassFlux(const label i, const surfaceScalarField& flux);
        
        //- Return the diffusion wall heat flux
        surfaceScalarField getDiffusiveWallHeatFlux() const;

        //- Read transportProperties dictionary
        virtual bool read() = 0;


        // Access

            //- Return the state of the solvingForX_ boolean
            inline bool solvingForX() const
            {
                return solvingForX_;
            }

            //- Return the state of the useNonCorrected_ boolean
            inline bool useNonCorrected() const
            {
                return useNonCorrected_;
            }

            //- Return the const binary diffusion coefficient field
            inline const volScalarField& Dij(label i, label j) const;

            //- Return the const product of the mixture density by the
            //  effective diffusion coefficient field for species i in a
            //  mixture
            inline const volScalarField& rhoD(label i) const;

            //- Return the const effective diffusion coefficient field for
            //  species i in a mixture
            inline const volScalarField& Di(label i) const;

            //- Return the table of species
            inline const speciesTable& species() const;

            //- Return the const mass flux field for species i
            inline const surfaceScalarField& speciesMassFlux
            (
                const label i
            ) const;

            //- Return the constant non-corrected diffusion mass flux field for
            //  a specie given by index
            inline const volVectorField& J(const label i) const;

            //- Return the constant diffusion mass flux field due to the
            //  gradient of pressure for a specie given by index
            inline const volVectorField& JGradp(const label i) const;

            //- Return the constant diffusion mass flux field due to the
            //  gradient of temperature for a specie given by index
            inline const volVectorField& JGradT(const label i) const;

            //- Return the const corrected diffusion mass flux field for a
            //  specie given by index
            inline const volVectorField Jcorr(const label i) const;

            //- Return the constant sum of the diffusion mass fluxes field
            inline const volVectorField& sumDiffusiveFluxes() const;

            //- Return the constant diffusion velocity field for a specie
            //  given by index
            inline const volVectorField& diffusionVelocity(const label i) const;

            //- Molar weight [kg/mol] for a specie given by index
            virtual scalar W(label i) const = 0;

            //- Absolute enthalpy [J/kg] for a specie given by index and
            //  temperature
            virtual scalar ha
            (
                label i,
                scalar p,
                scalar Tt,
                scalar Tv
            ) const = 0;
            
            //- Sensible enthalpy [J/kg] for a specie given by index and
            //  temperature
            virtual scalar hs
            (
                label i,
                scalar p,
                scalar Tt,
                scalar Tv
            ) const = 0;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "multiSpeciesTransportModelI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
