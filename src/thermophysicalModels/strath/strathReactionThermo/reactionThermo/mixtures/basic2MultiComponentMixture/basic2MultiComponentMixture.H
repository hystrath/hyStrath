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
    Foam::basic2MultiComponentMixture

Description
    Multi-component mixture. Provides a list of mass fraction fields and helper
    functions to query mixture composition.

SourceFiles
    basic2MultiComponentMixture.C

\*---------------------------------------------------------------------------*/

#ifndef basic2MultiComponentMixture_H
#define basic2MultiComponentMixture_H

#include "volFields.H"
#include "PtrList.H"
#include "speciesTable.H"

#include "DynamicList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class basic2MultiComponentMixture Declaration
\*---------------------------------------------------------------------------*/

class basic2MultiComponentMixture
{

protected:

    // Protected data

        //- Table of specie names
        speciesTable species_;
        
        //- Table of molecular specie names
        speciesTable molecules_;
        
        //- List of molecular specie indices
        labelList moleculeIds_;
        
        //- Table of atomic specie names
        speciesTable atoms_;
        
        //- List of atomic specie indices
        labelList atomIds_;
        
        //- Table of ionic specie names
        speciesTable ions_;
        
        //- List of ionic specie indices
        labelList ionIds_;
        
        //- Table of specie names, excluding free-electrons
        speciesTable heavySpecies_;
        
        //- List of specie indices, excluding free-electrons
        labelList heavySpeciesIds_;
        
        //- Position of the free-electrons in the list of species, if present
        label electronId_;

        //- List of specie active flags
        List<bool> active_;

        //- Table of vibrationally-excited specie names whose vib. eq. are
        //  being solved
        speciesTable solvedVibEqSpecies_;

        //- Table associating the vibro-electronic temperature of one specific
        //  molecule to atoms, ions, and electrons
        PtrList<label> vibTempAssociativity_;

        //- Species mass fractions
        PtrList<volScalarField> Y_;

        //- Species molar fractions
        PtrList<volScalarField> X_;

        //- Species number densities
        PtrList<volScalarField> nD_;

        //- Species number pressures
        PtrList<volScalarField> pP_;

        //- Species partial densities
        PtrList<volScalarField> pD_;

        //- Species vibrational (or vibro-electronic) temperatures
        PtrList<volScalarField> spTv_;

        //- Species vibrational (or vibro-electronic) temperatures
        //  for the different modes
        //PtrList<PtrList<volScalarField> > spmodeTv_; // TODO ABORTIVE WORK

        //- Species vibrational energies
        PtrList<volScalarField> hev_;

        //- Species electronic energies
        PtrList<volScalarField> heel_;

        //- Species vibro-electronic energies
        PtrList<volScalarField> hevel_;

        //- Species enthalpies
        PtrList<volScalarField> h_;

        //- Species vib-electronic energies
        //  for the different modes
        //PtrList<PtrList<volScalarField> > modehevel_; // TODO ABORTIVE WORK

        //- Energy field [J/kg]
        volScalarField* e_;

        //- Degrees of freedom associated with the rotational energy mode
        PtrList<volScalarField> zetaRot_;

        //- Degrees of freedom associated with the vibrational energy mode
        PtrList<volScalarField> zetaVib_;

        //- Degrees of freedom associated with the vibrational energy mode
        //  for the different modes
        //PtrList<PtrList<volScalarField> > modezetaVib_; // TODO ABORTIVE WORK

        //- Degrees of freedom associated with the vibrational energy mode
        PtrList<volScalarField> zetaElec_;

        //- Ouputs the desired fields
        bool writenD_, writepD_, writeX_, writepP_, writezeta_, writehev_;

        //- Mixture molecular weight [g/mol]
        volScalarField Wmix_;


    // Protected Member Functions

        //- Return the enthalpy/internal energy field boundary types
        //  by interrogating the temperature field boundary types
        wordList hev2BoundaryTypes
        (
            const label speciei,
            const volScalarField::Boundary& tbf1T,
            const bool downgradeToSingleTemperature = false
        );

        //- Return the enthalpy/internal energy field boundary base types
        //  by interrogating the temperature field boundary types
        wordList hev2BoundaryBaseTypes
        (
            const label speciei,
            const volScalarField::Boundary& tbf1T,
            const bool downgradeToSingleTemperature = false
        );

        //- Return the enthalpy/internal energy field boundary types
        //  by interrogating the temperature field boundary types
        wordList he2BoundaryTypes
        (
            const volScalarField::Boundary& tbf,
            const bool downgradeToSingleTemperature = false
        );

        //- Return the enthalpy/internal energy field boundary base types
        //  by interrogating the temperature field boundary types
        wordList he2BoundaryBaseTypes
        (
            const volScalarField::Boundary& tbf
        );


public:

    //- The base class of the mixture
    typedef basic2MultiComponentMixture basicMixtureType;

    // Constructors

        //- Construct from dictionary and mesh
        basic2MultiComponentMixture
        (
            const dictionary&,
            const wordList& specieNames,
            const fvMesh&
        );


    //- Destructor
    virtual ~basic2MultiComponentMixture()
    {
        delete e_;
    }


    // Member functions

        //- Return the table of species
        const speciesTable& species() const
        {
            return species_;
        }
        
        //- Return the table of molecular species
        const speciesTable& molecules() const
        {
            return molecules_;
        }
        
        //- Return the list of molecular species indices
        const labelList& moleculeIds() const
        {
            return moleculeIds_;
        }
        
        //- Return the moli^th molecule species index in species()
        const label& moleculeIds(const label moli) const
        {
            return moleculeIds_[moli];
        }
        
        //- Return the table of atomic species
        const speciesTable& atoms() const
        {
            return atoms_;
        }
        
        //- Return the list of atomic species indices
        const labelList& atomIds() const
        {
            return atomIds_;
        }
        
        //- Return the table of ionic species
        const speciesTable& ions() const
        {
            return ions_;
        }
        
        //- Return the list of ionic species indices
        const labelList& ionIds() const
        {
            return ionIds_;
        }
        
        //- Return the table of heavy species
        const speciesTable& heavySpecies() const
        {
            return heavySpecies_;
        }
        
        //- Return the list of heavy species indices
        const labelList& heavySpeciesIds() const
        {
            return heavySpeciesIds_;
        }
        
        //- Return the the speciei^th heavy species index in species()
        const label& heavySpeciesIds(const label speciei) const
        {
            return heavySpeciesIds_[speciei];
        }
        
        //- Return the position of the free-electrons in the list of species,
        //  if present
        const label& electronId() const
        {
            return electronId_;
        }

        //- Return the table of species
        const speciesTable& solvedVibEqSpecies() const
        {
            return solvedVibEqSpecies_;
        }

        //- Return the energy field
        inline volScalarField& e();

        //- Return the const energy field
        inline const volScalarField& e() const;

        //- Return the mass-fraction fields
        inline PtrList<volScalarField>& Y();

        //- Return the const mass-fraction fields
        inline const PtrList<volScalarField>& Y() const;

        //- Return the mass-fraction field for a specie given by index
        inline volScalarField& Y(const label i);

        //- Return the const mass-fraction field for a specie given by index
        inline const volScalarField& Y(const label i) const;

        //- Return the mass-fraction field for a specie given by name
        inline volScalarField& Y(const word& specieName);

        //- Return the const mass-fraction field for a specie given by name
        inline const volScalarField& Y(const word& specieName) const;

      //- X
        //- Return the molar-fraction fields
        inline PtrList<volScalarField>& X();

        //- Return the const molar-fraction fields
        inline const PtrList<volScalarField>& X() const;

        //- Return the molar-fraction field for a specie given by index
        inline volScalarField& X(const label i);

        //- Return the const molar-fraction field for a specie given by index
        inline const volScalarField& X(const label i) const;

        //- Return the molar-fraction field for a specie given by name
        inline volScalarField& X(const word& specieName);

        //- Return the const molar-fraction field for a specie given by name
        inline const volScalarField& X(const word& specieName) const;

      //- nD
        //- Return the number density fields
        inline PtrList<volScalarField>& nD();

        //- Return the const number density fields
        inline const PtrList<volScalarField>& nD() const;

        //- Return the number density field for a specie given by index
        inline volScalarField& nD(const label i);

        //- Return the const number density field for a specie given by index
        inline const volScalarField& nD(const label i) const;

        //- Return the number density field for a specie given by name
        inline volScalarField& nD(const word& specieName);

        //- Return the const number density field for a specie given by name
        inline const volScalarField& nD(const word& specieName) const;

      //- pP
        //- Return the partial pressures fields
        inline PtrList<volScalarField>& pP();

        //- Return the const partial pressures fields
        inline const PtrList<volScalarField>& pP() const;

        //- Return the partial pressures field for a specie given by index
        inline volScalarField& pP(const label i);

        //- Return the const partial pressures field for a specie given by index
        inline const volScalarField& pP(const label i) const;

        //- Return the partial pressures field for a specie given by name
        inline volScalarField& pP(const word& specieName);

        //- Return the const partial pressures field for a specie given by name
        inline const volScalarField& pP(const word& specieName) const;

      //- pD
        //- Return the partial densities fields
        inline PtrList<volScalarField>& pD();

        //- Return the const partial densities fields
        inline const PtrList<volScalarField>& pD() const;

        //- Return the partial densities field for a specie given by index
        inline volScalarField& pD(const label i);

        //- Return the const partial densities field for a specie given by index
        inline const volScalarField& pD(const label i) const;

        //- Return the partial densities field for a specie given by name
        inline volScalarField& pD(const word& specieName);

        //- Return the const partial densities field for a specie given by name
        inline const volScalarField& pD(const word& specieName) const;

      //- hev
        //- Return the vibrational energy fields
        inline PtrList<volScalarField>& hev();

        //- Return the const vibrational energy fields
        inline const PtrList<volScalarField>& hev() const;

        //- Return the vibrational energy field for a specie given by index
        inline volScalarField& hev(const label i);

        //- Return the const vibrational energy field for a specie given by
        //  index
        inline const volScalarField& hev(const label i) const;

        //- Return the vibrational energy field for a specie given by name
        inline volScalarField& hev(const word& specieName);

        //- Return the const vibrational energy field for a specie given by name
        inline const volScalarField& hev(const word& specieName) const;

      //- heel
        //- Return the electronic energy fields
        inline PtrList<volScalarField>& heel();

        //- Return the const electronic energy fields
        inline const PtrList<volScalarField>& heel() const;

        //- Return the electronic energy field for a specie given by index
        inline volScalarField& heel(const label i);

        //- Return the const electronic energy field for a specie given by index
        inline const volScalarField& heel(const label i) const;

        //- Return the electronic energy field for a specie given by name
        inline volScalarField& heel(const word& specieName);

        //- Return the const electronic energy field for a specie given by name
        inline const volScalarField& heel(const word& specieName) const;

      //- hevel
        //- Return the vibro-electronic energy fields
        inline PtrList<volScalarField>& hevel();

        //- Return the const vibro-electronic energy fields
        inline const PtrList<volScalarField>& hevel() const;

        //- Return the vibro-electronic energy field for a specie given by index
        inline volScalarField& hevel(const label i);

        //- Return the const vibro-electronic energy field for a specie given
        //  by index
        inline const volScalarField& hevel(const label i) const;

        //- Return the vibro-electronic energy field for a specie given by name
        inline volScalarField& hevel(const word& specieName);

        //- Return the const vibro-electronic energy field for a specie given
        //  by name
        inline const volScalarField& hevel(const word& specieName) const;

      //- h
        //- Return the enthalpy fields
        inline PtrList<volScalarField>& enthalpy();

        //- Return the const enthalpy fields
        inline const PtrList<volScalarField>& enthalpy() const;

        //- Return the enthalpy field for a specie given by index
        inline volScalarField& enthalpy(const label i);

        //- Return the const enthalpy field for a specie given by index
        inline const volScalarField& enthalpy(const label i) const;

        //- Return the enthalpy field for a specie given by name
        inline volScalarField& enthalpy(const word& specieName);

        //- Return the const enthalpy field for a specie given by name
        inline const volScalarField& enthalpy(const word& specieName) const;

      //- hevel for a specific mode m TODO ABORTIVE WORK
        /*//- Return the vib-electronic energy fields
        inline PtrList<PtrList<volScalarField> >& hevel_mode();

        //- Return the const vib-electronic energy fields
        inline const PtrList<PtrList<volScalarField> >& hevel_mode() const;

        //- Return the vib-electronic energy field for a specie given by index
        inline volScalarField& hevel_mode(const label i, const label m);

        //- Return the const vib-electronic energy field for a specie given by
        //  index
        inline const volScalarField& hevel_mode
        (
            const label i,
            const label m
        ) const;

        //- Return the electronic vib-energy field for a specie given by name
        inline volScalarField& hevel_mode
        (
            const word& specieName,
            const label i
        );

        //- Return the const electronic vib-energy field for a specie given by
        //  name
        inline const volScalarField& hevel_mode
        (
            const word& specieName,
            const label m
        ) const;*/

      //- Tv
        //- Return the vibrational temperature fields
        inline PtrList<volScalarField>& Tv();

        //- Return the const vibrational temperature fields
        inline const PtrList<volScalarField>& Tv() const;

        //- Return the vibrational temperature field for a specie given by index
        inline volScalarField& Tv(const label i);

        //- Return the const vibrational temperature field for a specie given
        //  by index
        inline const volScalarField& Tv(const label i) const;

        //- Return the vibrational temperature field for a specie given by name
        inline volScalarField& Tv(const word& specieName);

        //- Return the const vibrational temperature field for a specie given
        //  by name
        inline const volScalarField& Tv(const word& specieName) const;

      //- Tv for a specific mode m  TODO ABORTIVE WORK
        /*//- Return the vibrational temperature fields
        inline PtrList<PtrList<volScalarField> >& Tv_mode();

        //- Return the const vibrational temperature fields
        inline const PtrList<PtrList<volScalarField> >& Tv_mode() const;

        //- Return the vibrational temperature field for a specie given by index
        inline volScalarField& Tv_mode(const label i, const label m);

        //- Return the const vibrational temperature field for a specie given by
        //  index
        inline const volScalarField& Tv_mode
        (
            const label i,
            const label m
        ) const;

        //- Return the vibrational temperature field for a specie given by name
        inline volScalarField& Tv_mode(const word& specieName, const label m);

        //- Return the const vibrational temperature field for a specie given by
        //  name
        inline const volScalarField& Tv_mode
        (
            const word& specieName,
            const label m
        ) const;*/

      //- zetar
        //- Return the no of deg of freedom / rot energy mode fields
        inline PtrList<volScalarField>& zetar();

        //- Return the const no of deg of freedom / rot energy mode fields
        inline const PtrList<volScalarField>& zetar() const;

        //- Return the no of deg of freedom / rot energy mode field for a
        //  specie given by index
        inline volScalarField& zetar(const label i);

        //- Return the const no of deg of freedom / rot energy mode field for a
        //  specie given by index
        inline const volScalarField& zetar(const label i) const;

        //- Return the no of deg of freedom / rot energy mode field for a
        //  specie given by name
        inline volScalarField& zetar(const word& specieName);

        //- Return the const no of deg of freedom / rot energy mode field for a
        //  specie given by name
        inline const volScalarField& zetar(const word& specieName) const;

      //- zetav
        //- Return the no of deg of freedom / vib energy mode fields
        inline PtrList<volScalarField>& zetav();

        //- Return the const no of deg of freedom / vib energy mode fields
        inline const PtrList<volScalarField>& zetav() const;

        //- Return the no of deg of freedom / vib energy mode field for a
        //  specie given by index
        inline volScalarField& zetav(const label i);

        //- Return the const no of deg of freedom / vib energy mode field for a
        //  specie given by index
        inline const volScalarField& zetav(const label i) const;

        //- Return the no of deg of freedom / vib energy mode field for a
        //  specie given by name
        inline volScalarField& zetav(const word& specieName);

        //- Return the const no of deg of freedom / vib energy mode field for a
        //  specie given by name
        inline const volScalarField& zetav(const word& specieName) const;

      //- zetav  for a specific mode m TODO ABORTIVE WORK
        /*//- Return the no of deg of freedom / vib energy mode fields
        inline PtrList<PtrList<volScalarField> >& zetav_mode();

        //- Return the const no of deg of freedom / vib energy mode fields
        inline const PtrList<PtrList<volScalarField> >& zetav_mode() const;

        //- Return the no of deg of freedom / vib energy mode field for a
        //  specie given by index
        inline volScalarField& zetav_mode(const label i, const label m);

        //- Return the const no of deg of freedom / vib energy mode field for a
        //  specie given by index
        inline const volScalarField& zetav_mode
        (
            const label i,
            const label m
        ) const;

        //- Return the no of deg of freedom / vib energy mode field for a
        //  specie given by name
        inline volScalarField& zetav_mode
        (
            const word& specieName,
            const label m
        );

        //- Return the const no of deg of freedom / vib energy mode field for a
        //  specie given by name
        inline const volScalarField& zetav_mode
        (
            const word& specieName,
            const label m
        ) const;*/

      //- zetael
        //- Return the no of deg of freedom / elec energy mode fields
        inline PtrList<volScalarField>& zetael();

        //- Return the const no of deg of freedom / elec energy mode fields
        inline const PtrList<volScalarField>& zetael() const;

        //- Return the no of deg of freedom / elec energy mode field for a
        //  specie given by index
        inline volScalarField& zetael(const label i);

        //- Return the const no of deg of freedom / elec energy mode field for
        //  a specie given by index
        inline const volScalarField& zetael(const label i) const;

        //- Return the no of deg of freedom / elec energy mode field for a
        //  specie given by name
        inline volScalarField& zetael(const word& specieName);

        //- Return the const no of deg of freedom / elec energy mode field for
        //  a specie given by name
        inline const volScalarField& zetael(const word& specieName) const;

        //- Return the non-const mixture molecular weight field
        inline volScalarField& Wmix();

        //- Return the const mixture molecular weight field
        inline const volScalarField& Wmix() const;

        //- Does the mixture include this specie?
        inline bool contains(const word& specieName) const;
        
        //- Is this species a molecular specie of the mixture?
        inline bool isMolecule(const word& specieName) const;
        
        //- Does this species index point to a molecular specie of the mixture?
        inline bool isMolecule(const label speciei) const;
        
        //- Is this species an atomic specie of the mixture?
        inline bool isAtom(const word& specieName) const;
        
        //- Does this species index point to an atomic specie of the mixture?
        inline bool isAtom(const label speciei) const;
        
        //- Is this species an ionic specie of the mixture?
        inline bool isIon(const word& specieName) const;
        
        //- Does this species index point to an ionic specie of the mixture?
        inline bool isIon(const label speciei) const;
        
        //- Is this species a neutral specie of the mixture?
        inline bool isNeutral(const word& specieName) const;
        
        //- Does this species index point to a neutral specie of the mixture?
        inline bool isNeutral(const label speciei) const;

        //- Does the species have a vibration rate equation?
        inline bool containsVib(const word& specieName) const;

        //- Return true for active species
        inline bool active(label speciei) const;

        //- Return the bool list of active species
        inline const List<bool>& active() const;

        //- Set speciei active
        inline void setActive(label speciei);

        //- Set speciei inactive
        inline void setInactive(label speciei);

        inline scalar fres(const scalar ft, const scalar stoicRatio) const;

        inline tmp<volScalarField> fres
        (
            const volScalarField& ft,
            const dimensionedScalar& stoicRatio
        ) const;

        // TODO VINCENT: transform in inline fcts
        //- Mixture specific gas constant for cell-set [J/kg/K]
        virtual scalar RspecificMixture(const label celli) const = 0;

        //- Mixture specific gas constant for patch [J/kg/K]
        virtual scalar RspecificMixture
        (
            const label patchi,
            const label facei
        ) const = 0;

        //- Mixture molecular weight for cell-set [g/mol]
        virtual scalar molWeightMixture
        (
            const label celli,
            const bool fromMolarFrac = false
        ) const = 0;

        //- Mixture molecular weight for patch [g/mol]
        virtual scalar molWeightMixture
        (
            const label patchi,
            const label facei,
            const bool fromMolarFrac = false
        ) const = 0;

        //- Mixture molecular weight [g/mol]
        virtual tmp<volScalarField> molWeightMixture() const = 0;

        //- Update values of molar-fractions from mass-fractions for cell-set
        virtual scalar molarFraction
        (
            const label speciei,
            const scalar Yi,
            const label celli
        ) = 0;

        //- Update values of molar-fractions from mass-fractions for patch
        virtual scalar molarFraction
        (
            const label speciei,
            const scalar Yi,
            const label patchi,
            const label facei
        ) = 0;

        //- Update values of mass-fractions from molar-fractions during the
        //  initialisation
        virtual scalar massFractionFromMolarFraction
        (
            const label speciei,
            const scalar Xi,
            const scalar Wmix
        ) = 0;

        //- Update values of mass-fractions from partial densities
        virtual scalar massFractionFromPartialDensity
        (
            const scalar rhoi,
            const scalar p,
            const scalar Tt
        ) = 0;

        //- Update values of number densities from mass fractions
        virtual scalar numberDensity
        (
            const label speciei,
            const scalar Yi,
            const scalar rho
        ) = 0;

        //- Update values of partial pressures from molar fractions
        virtual scalar partialPressure(const scalar Xi, const scalar p) = 0;

        //- Update values of partial pressures from the equation of state
        virtual scalar partialPressureEoS
        (
            const label speciei,
            const scalar rhoi,
            const scalar Ti
        ) = 0;

        //- Update values of partial densities from mass fractions
        virtual scalar partialDensity(const scalar Yi, const scalar rho) = 0;


        // Per specie properties

            //- Number of moles []
            virtual scalar nMoles(const label specieI) const = 0;

            //- Molecular weight [kg/kmol]
            virtual scalar W(const label specieI) const = 0;

            //- Molecular weight [kg/kmol]
            virtual scalar W(const word& specieName) const = 0;

            //- Particle Type
            virtual label particleType(const label speciei) const = 0;

            //- Particle Charge
            virtual label particleCharge(const label speciei) const = 0;

            //- Diameter of the particle
            virtual scalar diameter(const label speciei) const = 0;

            //- Vibrational list with degeneracy and characteristic vibrational
            //  temperature [K] for each vib. mode
            virtual DynamicList<scalar> vibrationalList
            (
                const label speciei
            ) const = 0;

            //- Temperature-viscosity coefficient [-]
            virtual scalar omega(const label speciei) const = 0;

            //- Number of characteristic vibrational temperatures
            virtual label noVibrationalTemp(const label speciei) const = 0;

            //- Number of electronic energy levels
            virtual label noElectronicLevels(const label speciei) const = 0;

            //- Fraction of the first ionisation energy of a given species
            //  [J/kg]
            virtual scalar iHat(const label speciei) const = 0;

            //- Return the element of the table of vibro-electronic temperature
            //  associations for a specie given by index
            inline const label& vibTempAssociativity(const label speciei) const;

            //- Return the element of the table of vibro-electronic temperature
            //  associations for a specie given by name
            inline const label& vibTempAssociativity
            (
                const word& specieName
            ) const;


        // Per specie thermo properties

            //- Heat capacity at constant pressure [J/(kg K)]
            virtual scalar Cp
            (
                const label speciei,
                const scalar p,
                const scalar T
            ) const = 0;

            //- Heat capacity at constant volume [J/(kg K)]
            virtual scalar Cv
            (
                const label speciei,
                const scalar p,
                const scalar T
            ) const = 0;

            //- Absolute enthalpy [J/kg]
            virtual scalar Ha
            (
                const label speciei,
                const scalar p,
                const scalar T
            ) const = 0;

            //- Sensible enthalpy [J/kg]
            virtual scalar Hs
            (
                const label speciei,
                const scalar p,
                const scalar T
            ) const = 0;

            //- Chemical enthalpy [J/kg]
            virtual scalar Hc(const label speciei) const = 0;

            //- Entropy [J/(kg K)]
            virtual scalar S
            (
                const label speciei,
                const scalar p,
                const scalar T
            ) const = 0;

            //- Sensible internal energy [J/kg]
            virtual scalar Es
            (
                const label speciei,
                const scalar p,
                const scalar T
            ) const = 0;

            //- Gibbs free energy [J/kg]
            virtual scalar G
            (
                const label speciei,
                const scalar p,
                const scalar T
            ) const = 0;

            //- Helmholtz free energy [J/kg]
            virtual scalar A
            (
                const label speciei,
                const scalar p,
                const scalar T
            ) const = 0;

            //- Heat capacity at constant pressure [J/(kg K)]
            virtual scalar Cp
            (
                const label speciei,
                const scalar p,
                const scalar Tt,
                const scalar Tv
            ) const = 0;

            //- Heat capacity at constant volume [J/(kg K)]
            virtual scalar Cv
            (
                const label speciei,
                const scalar p,
                const scalar Tt,
                const scalar Tv
            ) const = 0;

            //- Heat capacity at constant volume [J/(kg K)]
            //  trans-rotational energy mode
            virtual scalar Cv_t
            (
                const label speciei,
                const scalar p,
                const scalar Tt
            ) const = 0;

            //- Heat capacity at constant volume [J/(kg K)]
            //  vibrational energy mode
            virtual scalar Cv_v
            (
                const label speciei,
                const scalar p,
                const scalar Tv
            ) const = 0;

            //- Heat capacity at constant volume [J/(kg K)]
            //  electronic energy mode
            virtual scalar Cv_el
            (
                const label speciei,
                const scalar p,
                const scalar Tv
            ) const = 0;

            //- Heat capacity at constant volume [J/(kg K)]
            //  vib-electronic energy mode
            virtual scalar Cv_vel
            (
                const label speciei,
                const scalar p,
                const scalar Tv
            ) const = 0;

            //- Heat capacity at constant pressure [J/(kg K)]
            //  trans-rotational energy mode
            virtual scalar Cp_t
            (
                const label speciei,
                const scalar p,
                const scalar Tt
            ) const = 0;

            //- Heat capacity at constant pressure [J/(kg K)]
            //  vibrational energy mode
            virtual scalar Cp_v
            (
                const label speciei,
                const scalar p,
                const scalar Tv
            ) const = 0;

            //- Heat capacity at constant pressure [J/(kg K)]
            //  electronic energy mode
            virtual scalar Cp_el
            (
                const label speciei,
                const scalar p,
                const scalar Tv
            ) const = 0;

            //- Heat capacity at constant pressure [J/(kg K)]
            //  vibro-electronic energy mode
            virtual scalar Cp_vel
            (
                const label speciei,
                const scalar p,
                const scalar Tv
            ) const = 0;

            //- Absolute enthalpy [J/kg]
            virtual scalar Ha
            (
                const label specieI,
                const scalar p,
                const scalar Tt,
                const scalar Tv
            ) const = 0;

            //- Sensible enthalpy [J/kg]
            virtual scalar Hs
            (
                const label specieI,
                const scalar p,
                const scalar Tt,
                const scalar Tv
            ) const = 0;

            //- Sensible enthalpy, trans-rotational [J/kg]
            virtual scalar Hts
            (
                const label speciei,
                const scalar p,
                const scalar Tt
            ) const = 0;

            //- Sensible enthalpy, vibrational [J/kg]
            virtual scalar Hvs
            (
                const label speciei,
                const scalar p,
                const scalar Tv
            ) const = 0;

            //- Sensible enthalpy, electronic [J/kg]
            virtual scalar Hels
            (
                const label speciei,
                const scalar p,
                const scalar Tv
            ) const = 0;

            //- Sensible enthalpy, vibro-electronic [J/kg]
            virtual scalar Hvels
            (
                const label speciei,
                const scalar p,
                const scalar Tv
            ) const = 0;

            //- Entropy [J/(kg K)]
            virtual scalar S
            (
                const label specieI,
                const scalar p,
                const scalar Tt,
                const scalar Tv
            ) const = 0;

            //- Sensible internal energy [J/kg]
            virtual scalar Es
            (
                const label specieI,
                const scalar p,
                const scalar Tt,
                const scalar Tv
            ) const = 0;

            //- Sensible internal energy, trans-rotational [J/kg]
            virtual scalar Ets
            (
                const label speciei,
                const scalar p,
                const scalar Tt
            ) const = 0;

            //- Sensible internal energy, vibrational [J/kg]
            virtual scalar Evs
            (
                const label speciei,
                const scalar p,
                const scalar Tv
            ) const = 0;

            //- Sensible internal energy, electronic [J/kg]
            virtual scalar Eels
            (
                const label speciei,
                const scalar p,
                const scalar Tv
            ) const = 0;

            //- Sensible internal energy, vib-electronic [J/kg]
            virtual scalar Evels
            (
                const label speciei,
                const scalar p,
                const scalar Tv
            ) const = 0;

            //- Trans-rotational Enthalpy/Internal energy [J/kg]
            virtual scalar HEt
            (
                const label speciei,
                const scalar p,
                const scalar Tt
            ) const = 0;

            //- Vibrational Enthalpy/Internal energy [J/kg]
            virtual scalar HEv
            (
                const label speciei,
                const scalar p,
                const scalar Tv
            ) const = 0;

            //- Electronic Enthalpy/Internal energy [J/kg]
            virtual scalar HEel
            (
                const label speciei,
                const scalar p,
                const scalar Tv
            ) const = 0;

            //- Vibro-electronic Enthalpy/Internal energy [J/kg]
            virtual scalar HEvel
            (
                const label speciei,
                const scalar p,
                const scalar Tv
            ) const = 0;

            //- Vibro-electronic Enthalpy/Internal energy [J/kg]
            //  for the vibrational mode m TODO ABORTIVE WORK
            /*virtual scalar HEvel_mode
            (
                const label speciei,
                const label mode,
                const scalar p,
                const scalar Tv
            ) const = 0;*/

            //- Heat capacity at constant volume for patch[J/kg/K]
            //  vibro-electronic energy mode
            virtual tmp<scalarField> Cv_vel
            (
                const label speciei,
                const scalarField& p,
                const scalarField& Tv,
                const label patchi
            ) const = 0;

            //- Heat capacity at constant volume for patch[J/kg/K]
            //  vibro-electronic energy mode
            //  Interface with the energy boundary condition
            virtual tmp<scalarField> Cv_vel
            (
                const word specieName,
                const scalarField& p,
                const scalarField& Tv,
                const label patchi
            ) const
            {
                return Cv_vel
                (
                    species_[specieName],
                    p,
                    Tv,
                    patchi
                );
            }

            //- Vibro-electronic Enthalpy/Internal energy for cell-set [J/kg]
            virtual tmp<scalarField> hevel
            (
                const label speciei,
                const scalarField& p,
                const scalarField& Tv,
                const labelList& cells
            ) const = 0;

            //- Vibro-electronic Enthalpy/Internal energy for patch [J/kg]
            //  Interface with the energy boundary condition
            virtual tmp<scalarField> hevel
            (
                const word specieName,
                const scalarField& p,
                const scalarField& Tv,
                const labelList& cells
            ) const
            {
                return hevel
                (
                    species_[specieName],
                    p,
                    Tv,
                    cells
                );
            }

            //- Vibro-electronic Enthalpy/Internal energy for patch [J/kg]
            virtual tmp<scalarField> hevel
            (
                const label speciei,
                const scalarField& p,
                const scalarField& Tv,
                const label patchi
            ) const = 0;

            //- Vibro-electronic Enthalpy/Internal energy for patch [J/kg]
            //  Interface with the energy boundary condition
            virtual tmp<scalarField> hevel
            (
                const word specieName,
                const scalarField& p,
                const scalarField& Tv,
                const label patchi
            ) const
            {
                return hevel
                (
                    species_[specieName],
                    p,
                    Tv,
                    patchi
                );
            }

            //- Trans-rotational temperature from sensible internal energy
            //  given an initial T0
            virtual scalar TtHEt
            (
                const label speciei,
                const scalar het,
                const scalar p,
                const scalar Tt
            ) const = 0;

            //- Vibrational temperature from sensible internal energy
            //  given an initial T0
            virtual scalar TvHEv
            (
                const label speciei,
                const scalar hev,
                const scalar p,
                const scalar Tv
            ) const = 0;

            //- Vibro-electronic temperature from sensible internal energy
            //  given an initial T0
            virtual scalar TvelHEvel
            (
                const label speciei,
                const scalar hevel,
                const scalar p,
                const scalar Tvel
            ) const = 0;

            //- Vib-electronic temperature from sensible internal energy given an
            //  initial T0, for the vibrational mode m TODO ABORTIVE WORK
            /*virtual scalar TvelHEvel_mode
            (
                const label speciei,
                const label mode,
                const scalar hevel,
                const scalar p,
                const scalar Tvel
            ) const = 0;*/

            //- Degrees of freedom associated with the rotational mode [-]
            virtual scalar zetar
            (
                const label speciei,
                const scalar p,
                const scalar Tt,
                const scalar Tv
            ) const = 0;

            //- Degrees of freedom associated with the vibrational mode [-]
            virtual scalar zetav
            (
                const label speciei,
                const scalar p,
                const scalar Tv
            ) const = 0;

            //- Degrees of freedom associated with the vibrational mode [-]
            //  for the vibrational mode m TODO ABORTIVE WORK
            /*virtual scalar zetav_mode
            (
                const label speciei,
                const label mode,
                const scalar p,
                const scalar Tv
            ) const = 0;*/

            //- Degrees of freedom associated with the electronic mode [-]
            virtual scalar zetael
            (
                const label speciei,
                const scalar p,
                const scalar Tv
            ) const = 0;

            //- Gibbs free energy [J/kg]
            virtual scalar G
            (
                const label specieI,
                const scalar p,
                const scalar Tt,
                const scalar Tv
            ) const = 0;

            //- Helmholtz free energy [J/kg]
            virtual scalar A
            (
                const label specieI,
                const scalar p,
                const scalar Tt,
                const scalar Tv
            ) const = 0;

        // Per specie transport properties

            //- Compressibility rho/p [s^2/m^2]
            virtual scalar psi
            (
                const label speciei,
                const scalar p,
                const scalar T
            ) const = 0;

            //- Dynamic viscosity [kg/m/s]
            virtual scalar mu
            (
                const label specieI,
                const scalar p,
                const scalar T
            ) const = 0;

            //- Thermal conductivity [W/m/K]
            //  trans-rotational energy mode
            virtual scalar kappatr
            (
                const label speciei,
                const scalar p,
                const scalar Tt
            ) const = 0;

            //- Thermal conductivity [W/m/K]
            //  vibro-electronic energy mode
            virtual scalar kappave
            (
                const label speciei,
                const scalar p,
                const scalar Tt,
                const scalar Tve = 0
            ) const = 0;

            //- Thermal diffusivity [kg/m/s]
            //  trans-rotational energy mode
            virtual scalar alphatr
            (
                const label speciei,
                const scalar p,
                const scalar Tt
            ) const = 0;

            //- Thermal diffusivity [kg/m/s]
            //  vibro-electronic energy mode
            virtual scalar alphave
            (
                const label speciei,
                const scalar p,
                const scalar Tt,
                const scalar Tve = 0
            ) const = 0;

            //- Density [kg/m3]
            virtual scalar rho
            (
                const label specieI,
                const scalar p,
                const scalar T
            ) const = 0;


        //- Outputs the desired fields
        virtual void write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "basic2MultiComponentMixtureI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
