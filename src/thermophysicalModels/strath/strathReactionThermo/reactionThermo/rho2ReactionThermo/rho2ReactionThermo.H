/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

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
    Foam::rho2ReactionThermo

Description
    Foam::rho2ReactionThermo

SourceFiles
    rho2ReactionThermo.C

\*---------------------------------------------------------------------------*/

#ifndef rho2ReactionThermo_H
#define rho2ReactionThermo_H

#include "rho2Thermo.H"
#include "basic2MultiComponentMixture.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"

#include "constants.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class rho2ReactionThermo Declaration
\*---------------------------------------------------------------------------*/

class rho2ReactionThermo
:
    public rho2Thermo
{

protected:

    // Protected data
      
        //- Partial name of the thermo
        word partialThermoName_;
        
        //- Temperature limits of applicability of functions (for all species)
        scalar Tlow_, Thigh_;
        
        //- Minimum and maximum temperatures tolerated at patches
        scalar TlowPatches_, ThighPatches_;
        
        //- Maximum velocity tolerated at patches (magnitude for each component)
        scalar UhighPatches_;
        

    // Private Member Functions
        
        //- Translates the transport model name into the typedef name of the partial thermo
        word transportToTypedef(const word transportModel);
        
        //- Determines the temperature from an energy field F, its derivate dF, and the mass fractions YList
        inline scalar Temperature
        (
            scalar f,
            scalar p,
            scalar T0,
            const PtrList<scalar>& YList,
            scalar (Foam::basic2MultiComponentMixture::*F)(const label, const scalar, const scalar) const,
            scalar (Foam::basic2MultiComponentMixture::*dFdT)(const label, const scalar, const scalar) const,
            scalar (Foam::rho2ReactionThermo::*limit)(const scalar) const
        ) const;
        
        //- Determines the temperature from an energy field F, its derivate dF for a given specie i
        inline scalar Temperature
        (
            label i,
            scalar f,
            scalar p,
            scalar T0,
            scalar (Foam::basic2MultiComponentMixture::*F)(const label, const scalar, const scalar) const,
            scalar (Foam::basic2MultiComponentMixture::*dFdT)(const label, const scalar, const scalar) const,
            scalar (Foam::rho2ReactionThermo::*limit)(const scalar) const
        ) const;
        
        //- Overall temperature calculation for the single-temperature model
        inline scalar TEs // NEW VINCENT 16/06/2016
        (
            const scalar e,
            const scalar p,
            const scalar T0,
            const PtrList<scalar>& YList
        ) const;
        
        //- Trans-rotational temperature calculation for the two-temperature model
        inline scalar TtEts // NEW VINCENT 17/06/2016
        (
            const scalar e,
            const scalar p,
            const scalar T0,
            const PtrList<scalar>& YList
        ) const;
        
        //- Single-vibrational temperature calculation for the two-temperature model
        inline scalar TvelEvels // NEW VINCENT 17/06/2016
        (
            const scalar e,
            const scalar p,
            const scalar T0,
            const PtrList<scalar>& YList
        ) const;
        
        //- Multi-vibrational temperature calculation for the two-temperature model
        inline scalar TvelEvels // NEW VINCENT 19/08/2016
        (
            const label i,
            const scalar e,
            const scalar p,
            const scalar T0
        ) const;
        
        //- Limit the temperature to be in the range Tlow_ to Thigh_
        inline scalar limit(const scalar T) const;
        
        //- Calculate the chemistry fractions
        void correctChemFractions();
        
        //- Calculate the fields from DSMC inputs
        void calculateFromDSMC
        (
            PtrList<volScalarField>& fixedY,
            volScalarField& fixedTt,
            PtrList<volScalarField>& fixedTv,
            volScalarField& fixedRho,
            word& regionName
        );
        
        //- Calculate the thermo variables
        void calculate();
        
        //- Minimal calculation of the thermo variables
        void calculateLight();
        
        //- Calculate hevk energies of non-molecule particles only
        void calculateNonMoleculeHEVK();
        
        //- Initialise the fields
        void initialise();
        
        //- Minimal initialisation of the fields
        void initialiseLight();

public:

    //- Public constant
     
        //- Vibrational cut-off temperature
        static const scalar vibrationalCutOffTemp;
        
        //- Is there at least one cell for which the temperature is not
        //  within [Tlow, Thigh]?
        static bool temperatureFieldOutOfRange;
        
    //- Runtime type information
    TypeName("rho2ReactionThermo");


    //- Declare run-time constructor selection tables
    declareRunTimeSelectionTable
    (
        autoPtr,
        rho2ReactionThermo,
        fvMesh,
        (const fvMesh& mesh, const word& phaseName),
        (mesh, phaseName)
    );


    // Constructors

        //- Construct from mesh and phase name
        rho2ReactionThermo
        (
            const fvMesh&,
            const word& phaseName
        );


    // Selectors

        //- Standard selection based on fvMesh
        static autoPtr<rho2ReactionThermo> New
        (
            const fvMesh&,
            const word& phaseName=word::null
        );


    //- Destructor
    virtual ~rho2ReactionThermo();


    // Member functions

        //- Constant access to the ptach field values for energy boundary calculations
        static const rho2ReactionThermo& lookup2ReactionThermo(const fvPatchScalarField& pf);
        
        //- Return the partial name of the thermo package that has been loaded
        word partialThermoName();
        
        //- Return the composition of the multi-component mixture
        virtual basic2MultiComponentMixture& composition() = 0;

        //- Return the composition of the multi-component mixture
        virtual const basic2MultiComponentMixture& composition() const = 0;

        //- Update the fields from DSMC inputs
        virtual void correctFromDSMC
        (
            PtrList<volScalarField>& fixedY,
            volScalarField& fixedTt,
            PtrList<volScalarField>& fixedTv,
            volScalarField& fixedRho,
            word& regionName
        );
        
        //- Update properties
        virtual void correct2T();
        
        //- Minimal properties update
        virtual void correct2T_Light();
        
        //- Update hevk energies of non-molecule particles only
        virtual void correctNonMoleculeHEVK();
        
        //- Limit the velocity field at wall patches
        virtual void limitVelocityAtWallBoundary(volVectorField::GeometricBoundaryField& Ubdry);
        
        //- Update chemistry fractions only
        virtual void correctFractions();
        
        //- Update the overall temperature
        virtual void correctOverallTemperature();
        
        //- Initialise the 2T model
        virtual void initialise2T();
        
        //- Initialise the 2T model (light)
        virtual void initialise2T_Light();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "rho2ReactionThermoI.H"

#endif

// ************************************************************************* //
