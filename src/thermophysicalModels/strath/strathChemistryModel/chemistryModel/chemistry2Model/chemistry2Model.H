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
    Foam::chemistry2Model

Description
    Extends base chemistry model by adding a thermo package, and ODE functions.
    Introduces chemistry equation system and evaluation of chemical source
    terms.

SourceFiles
    chemistry2ModelI.H
    chemistry2Model.C

\*---------------------------------------------------------------------------*/

#ifndef chemistry2Model_H
#define chemistry2Model_H

#include "Reaction2.H"
#include "ODESystem.H"
#include "volFieldsFwd.H"
#include "simpleMatrix.H"
#include "DimensionedField.H"
#include "Switch.H"

#include "constants.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class fvMesh;

/*---------------------------------------------------------------------------*\
                      Class chemistry2Model Declaration
\*---------------------------------------------------------------------------*/

template<class CompType, class ThermoType>
class chemistry2Model
:
    public CompType,
    public ODESystem
{
    // Private Member Functions

        //- Disallow copy constructor
        chemistry2Model(const chemistry2Model&);

        //- Disallow default bitwise assignment
        void operator=(const chemistry2Model&);

        //- Solve the reaction system for the given time step
        //  of given type and return the characteristic time
        template<class DeltaTType>
        scalar solve(const DeltaTType& deltaT);

        //- Return the modified temperature when bounding (fixed min value) is
        //  required. Useful when equilibrium constant is a very small value
        //  cf Scalabrin's PhD thesis (2007), Eq. 3.13, p37
        //  from experience, needed for electron impact reactions
        scalar modTemperature(const scalar T) const;


protected:

    typedef ThermoType thermoType;

    // Private data

        //- Reference to the field of specie mass fractions
        PtrList<volScalarField>& Y_;

        //- Reactions
        const PtrList<Reaction2<ThermoType> >& reactions_;

        //- Thermodynamic data of the species
        const PtrList<ThermoType>& specieThermo_;

        //- Number of species
        label nSpecie_;

        //- Number of reactions
        label nReaction_;

        //- Temperature below which the reaction rates are assumed 0
        scalar Treact_;

        //- List of reaction rate per specie [kg/m3/s]
        PtrList<DimensionedField<scalar, volMesh> > RR_;

        //- List of forward reaction rate per specie [kg/m3/s]
        PtrList<DimensionedField<scalar, volMesh> > RRf_;

        //- Flag to indicate whether any of the ionised atoms N+ and O+ species
        //  are present in the mixture
        Switch containsIonisedAtoms_;
        
        //- Position of the ionised atoms N+ and O+ in the species list,
        //  respectively, if present
        FixedList<label, 2> posIonisedAtoms_;
        
        //- Position of the neutral atoms N and O in the species list,
        //  respectively, if needed (that is, if there are ionised atoms in the
        //  mixture)
        FixedList<label, 2> posNeutralAtoms_;

        // Park TTv model
        
          //- Factor in the preferential model, default value: 0.3
          PtrList<FixedList<scalar, 2> > preferentialFactor_;

        // CVDV model
        
          //- Reciprocal U factor
          scalar reciprocalUFactor_;

          //- Number of vibrational energy levels considered in the simple
          //  harmonic oscillator model
          PtrList<label> simpleHarmonicOscillatorVibCutOff_;


    // Protected Member Functions

        //- Write access to chemical source terms
        //  (e.g. for multi-chemistry model)
        inline PtrList<DimensionedField<scalar, volMesh> >& RR();


public:

    //- Runtime type information
    TypeName("chemistry2Model");


    // Constructors

        //- Construct from mesh
        chemistry2Model(const fvMesh& mesh);


    //- Destructor
    virtual ~chemistry2Model();


    // Member Functions

        //- The reactions
        inline const PtrList<Reaction2<ThermoType> >& reactions() const;

        //- Thermodynamic data of the species
        inline const PtrList<ThermoType>& specieThermo() const;

        //- The number of species
        inline label nSpecie() const;

        //- The number of reactions
        inline label nReaction() const;

        //- Temperature below which the reaction rates are assumed 0
        inline scalar Treact() const;

        //- Temperature below which the reaction rates are assumed 0
        inline scalar& Treact();

        //- dc/dt = omega, rate of change in concentration, for each species
        virtual tmp<scalarField> omega
        (
            const scalarField& c,
            const scalar T,
            const scalar Tv,
            const scalar p
        ) const;

        //- Return the reaction rate for reaction r and the reference
        //  species and charateristic times
        virtual scalar omega
        (
            const Reaction2<ThermoType>& R,
            const scalarField& c,
            const scalar T,
            const scalar Tv,
            const scalar p,
            scalar& pf,
            scalar& cf,
            label& lRef,
            scalar& pr,
            scalar& cr,
            label& rRef
        ) const;


        //- Return the reaction rate for reaction r and the reference
        //  species and charateristic times, overloaded for Tv[speciei]
        virtual scalar omega
        (
            const Reaction2<ThermoType>& r,
            const scalarField& c,
            const scalar T,
            const scalarList& spTv,
            const scalar p,
            scalar& pf,
            scalar& cf,
            label& lRef,
            scalar& pr,
            scalar& cr,
            label& rRef
        ) const;


        //- Return the reaction rate for iReaction and the reference
        //  species and charateristic times
        virtual scalar omegaI
        (
            label iReaction,
            const scalarField& c,
            const scalar T,
            const scalarList& spTv,
            const scalar p,
            scalar& pf,
            scalar& cf,
            label& lRef,
            scalar& pr,
            scalar& cr,
            label& rRef
        ) const;

        //- Calculates the reaction rates
        virtual void calculate();


        // Chemistry model functions (overriding abstract functions in
        // basic2ChemistryModel.H)

            //- Return const access to the chemical source terms for specie, i
            inline const DimensionedField<scalar, volMesh>& RR
            (
                const label i
            ) const;

            //- Return non const access to chemical source terms [kg/m/s3]
            virtual DimensionedField<scalar, volMesh>& RR
            (
                const label i
            );

            //- Return const access to the forward chemical source terms for
            //  specie, i
            inline const DimensionedField<scalar, volMesh>& RRf
            (
                const label i
            ) const;

            //- Return non const access to forward chemical source terms
            //  units: [kg/m/s3]
            virtual DimensionedField<scalar, volMesh>& RRf
            (
                const label i
            );

            //- Return reaction rate of the specieI in reactionI
            virtual tmp<DimensionedField<scalar, volMesh> > calculateRR
            (
                const label reactionI,
                const label specieI
            ) const;

            //- Solve the reaction system for the given time step
            //  and return the characteristic time
            virtual scalar solve(const scalar deltaT);

            //- Solve the reaction system for the given time step
            //  and return the characteristic time
            virtual scalar solve(const scalarField& deltaT);

            //- Return the chemical time scale
            virtual tmp<volScalarField> tc() const;

            //- Return source for enthalpy equation [kg/m/s3]
            virtual tmp<volScalarField> Sh() const;

            //- Return chemistry-vibrational coupling source term
            //  in the vibrational equation [kg/m/s3]
            virtual tmp<volScalarField> Scv() const;
            
            //- Return chemistry-vibrational coupling source term
            //  in the vibrational equation for specie i [kg/m/s3]
            virtual tmp<volScalarField> Scv(const label i) const;

            //- Return energy removed from free-electrons during impact
            //  ionisation reactions [kg/m/s3]
            virtual tmp<volScalarField> Seiir() const;
            
            //- Return energy removed from free-electrons during impact
            //  ionisation reactions for specie i [kg/m/s3]
            virtual tmp<volScalarField> Seiir(const label i) const;

            //- Return the heat release, i.e. enthalpy/sec [kg/m2/s3]
            virtual tmp<volScalarField> dQ() const;


        // ODE functions (overriding abstract functions in ODE2.H)

            //- Number of ODE's to solve
            virtual label nEqns() const;

            virtual void derivatives
            (
                const scalar t,
                const scalarField& c,
                scalarField& dcdt
            ) const;

            virtual void jacobian
            (
                const scalar t,
                const scalarField& c,
                scalarField& dcdt,
                scalarSquareMatrix& dfdc
            ) const;

            virtual void solve
            (
                scalarField &c,
                scalar& T,
                scalarList& spTv,
                scalar& p,
                scalar& deltaT,
                scalar& subDeltaT
            ) const;

            virtual void solve
            (
                scalarField& c,
                scalarField& cwfd,
                scalar& T,
                scalarList& spTv,
                scalar& p,
                scalar& deltaT,
                scalar& subDeltaT
            ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "chemistry2ModelI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "chemistry2Model.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
