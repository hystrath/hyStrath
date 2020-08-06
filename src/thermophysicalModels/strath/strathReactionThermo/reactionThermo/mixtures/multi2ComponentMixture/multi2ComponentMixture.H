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
    Foam::multi2ComponentMixture

Description
    Foam::multi2ComponentMixture

SourceFiles
    multi2ComponentMixture.C

\*---------------------------------------------------------------------------*/

#ifndef multi2ComponentMixture_H
#define multi2ComponentMixture_H

#include "basic2MultiComponentMixture.H"
#include "HashPtrTable.H"

#include "fundamentalConstants.H"
#include "constants.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class multi2ComponentMixture Declaration
\*---------------------------------------------------------------------------*/

template<class ThermoType>
class multi2ComponentMixture
:
    public basic2MultiComponentMixture
{
    // Private data

        //- Species data
        PtrList<ThermoType> speciesData_;

        //- Temporary storage for the cell/face mixture thermo data
        mutable ThermoType mixture_;

        //- Temporary storage for the volume weighted
        //  cell/face mixture thermo data
        mutable ThermoType mixtureVol_;


    // Private Member Functions

        //- Construct the species data from the given dictionary and return the
        //  data for the first specie to initialise the mixture thermo data
        const ThermoType& constructSpeciesData(const dictionary& thermoDict);

        //- Correct the mass fractions to sum to 1
        void correctMassFractions();

        //- Correct the associativity vibrational temperature table
        void correctVibTempAssociativity(const dictionary& dict);

        //- Fills the solvedVibEqSpecies table from the list of species and the
        //  vibTempAssociativity list
        void fillSolvedVibEqSpeciesTable();

        //- Fills a modified solvedVibEqSpecies table from the list of species
        //  that will later serves in the code
        void fillHackOfSolvedVibEqSpeciesTable();

        //- Construct as copy (not implemented)
        multi2ComponentMixture(const multi2ComponentMixture<ThermoType>&);


public:

    //- The type of thermodynamics this mixture is instantiated for
    typedef ThermoType thermoType;


    // Constructors

        //- Construct from dictionary, specie names, thermo database and mesh
        multi2ComponentMixture
        (
            const dictionary&,
            const wordList& specieNames,
            const HashPtrTable<ThermoType>& thermoData,
            const fvMesh&
        );

        //- Construct from dictionary and mesh
        multi2ComponentMixture(const dictionary&, const fvMesh&);


    //- Destructor
    virtual ~multi2ComponentMixture()
    {}


    // Member functions

        const ThermoType& cellMixture(const label celli) const;

        const ThermoType& patchFaceMixture
        (
            const label patchi,
            const label facei
        ) const;

        scalar cellMixture_HEt
        (
            const label celli,
            const scalar p,
            const scalar Tt
        ) const;

        scalar cellMixture_HEv
        (
            const label celli,
            const scalar p,
            const scalar Tv
        ) const;

        scalar cellMixture_HEel
        (
            const label celli,
            const scalar p,
            const scalar Tv
        ) const;

        scalar cellMixture_HEvel
        (
            const label celli,
            const scalar p,
            const scalar Tv
        ) const;

        scalar cellMixture_Cv_t
        (
            const label celli,
            const scalar p,
            const scalar Tt
        ) const;

        scalar cellMixture_Cp_t
        (
            const label celli,
            const scalar p,
            const scalar Tt
        ) const;

        scalar cellMixture_gamma
        (
            const label celli,
            const scalar p,
            const scalar Tt,
            const scalar Tv
        ) const;

        scalar cellMixture_kappa_tr
        (
            const label celli,
            const scalar p,
            const scalar Tt
        ) const;

        scalar cellMixture_kappa_ve
        (
            const label celli,
            const scalar p,
            const scalar Tt,
            const scalar Tv
        ) const;

        scalar cellMixture_alpha_tr
        (
            const label celli,
            const scalar p,
            const scalar Tt
        ) const;

        scalar cellMixture_alpha_ve
        (
            const label celli,
            const scalar p,
            const scalar Tt,
            const scalar Tv
        ) const;

        scalar patchFaceMixture_HEt
        (
            const label patchi,
            const label facei,
            const scalar p,
            const scalar Tt
        ) const;

        scalar patchFaceMixture_HEv
        (
            const label patchi,
            const label facei,
            const scalar p,
            const scalar Tv
        ) const;

        scalar patchFaceMixture_HEel
        (
            const label patchi,
            const label facei,
            const scalar p,
            const scalar Tv
        ) const;

        scalar patchFaceMixture_HEvel
        (
            const label patchi,
            const label facei,
            const scalar p,
            const scalar Tv
        ) const;

        scalar patchFaceMixture_Cv_t
        (
            const label patchi,
            const label facei,
            const scalar p,
            const scalar Tt
        ) const;

        scalar patchFaceMixture_Cp_t
        (
            const label patchi,
            const label facei,
            const scalar p,
            const scalar Tt
        ) const;

        scalar patchFaceMixture_gamma
        (
            const label patchi,
            const label facei,
            const scalar p,
            const scalar Tt,
            const scalar Tv
        ) const;

        scalar patchFaceMixture_kappa_tr
        (
            const label patchi,
            const label facei,
            const scalar p,
            const scalar Tt
        ) const;

        scalar patchFaceMixture_kappa_ve
        (
            const label patchi,
            const label facei,
            const scalar p,
            const scalar Tt,
            const scalar Tv
        ) const;

        scalar patchFaceMixture_alpha_tr
        (
            const label patchi,
            const label facei,
            const scalar p,
            const scalar Tt
        ) const;

        scalar patchFaceMixture_alpha_ve
        (
            const label patchi,
            const label facei,
            const scalar p,
            const scalar Tt,
            const scalar Tv
        ) const;

        //- Attempt to return the value of any function F for patch (not used)
        scalar patchFaceMixture
        (
            scalar (ThermoType::*F)(const scalar, const scalar) const,
            const label patchi,
            const label facei,
            const scalar p,
            const scalar T
        ) const;

        const ThermoType& cellVolMixture
        (
            const scalar p,
            const scalar T,
            const label celli
        ) const;

        const ThermoType& patchFaceVolMixture
        (
            const scalar p,
            const scalar T,
            const label patchi,
            const label facei
        ) const;

        //- Mixture molecular weight for cell-set [g/mol]
        virtual scalar molWeightMixture
        (
            const label celli,
            const bool fromMolarFrac = false
        ) const;

        //- Mixture molecular weight for patch [g/mol]
        virtual scalar molWeightMixture
        (
            const label patchi,
            const label facei,
            const bool fromMolarFrac = false
        ) const;
        
        //- Mixture molecular weight [g/mol]
        virtual tmp<volScalarField> molWeightMixture() const;

        //- Mixture specific gas constant for cell-set [J/kg/K]
        virtual scalar RspecificMixture(const label celli) const;

        //- Mixture specific gas constant for cell-set for patch [J/kg/K]
        virtual scalar RspecificMixture
        (
            const label patchi,
            const label facei
        ) const;

        //- Update values of molar-fractions from mass-fractions for cell-set
        virtual scalar molarFraction
        (
            const label speciei,
            const scalar Yi,
            const label celli
        );

        //- Update values of molar-fractions from mass-fractions for patch
        virtual scalar molarFraction
        (
            const label speciei,
            const scalar Yi,
            const label patchi,
            const label facei
        );

        //- Update values of mass-fractions from molar-fractions during the
        //  initialisation
        virtual scalar massFractionFromMolarFraction
        (
            const label speciei,
            const scalar Xi,
            const scalar Wmix
        );

        //- Update values of mass-fractions from partial densities
        virtual scalar massFractionFromPartialDensity
        (
            const scalar rhoi,
            const scalar p,
            const scalar Tt
        );

        //- Update values of number densities from mass fractions
        virtual scalar numberDensity
        (
            const label speciei,
            const scalar Yi,
            const scalar rho
        );

        //- Update values of partial pressures from molar fractions
        virtual scalar partialPressure(const scalar Xi, const scalar p);

        //- Update values of partial pressures from the equation of state
        virtual scalar partialPressureEoS
        (
            const label speciei,
            const scalar rhoi,
            const scalar Ti
        );

        //- Update values of partial densities from mass fractions
        virtual scalar partialDensity(const scalar Yi, const scalar rho);

        //- Return the raw specie thermodynamic data
        const PtrList<ThermoType>& speciesData() const
        {
            return speciesData_;
        }

        //- Read dictionary
        void read(const dictionary&);

        //- Return thermo based on index
        inline const ThermoType& getLocalThermo(const label speciei) const
        {
            return speciesData_[speciei];
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "multi2ComponentMixture.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
