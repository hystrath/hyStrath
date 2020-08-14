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
    Foam::multi2Thermo

Description
    Fundamental multi-temperature thermodynamic properties

SourceFiles
    multi2Thermo.C

\*---------------------------------------------------------------------------*/

#ifndef multi2Thermo_H
#define multi2Thermo_H

#include "fluid2Thermo.H" // NEW VINCENT 19/04/2016

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class multi2Thermo Declaration
\*---------------------------------------------------------------------------*/

class multi2Thermo
:
    public fluid2Thermo
{


protected:

    // Protected data

        //- Two-temperature dictionary
        IOdictionary twoTemperatureDictionary_;

        //- Overall temperature [K]
        volScalarField Tov_;

        //- Vibro-electronic temperature of the mixture [K]
        volScalarField Tv_;

        //- Energy field associated with the Overall mode [J/kg]
        volScalarField het_;

        //- Energy field associated with the vibrational mode [J/kg]
        volScalarField hevMix_;

        //- Energy field associated with the electronic mode [J/kg]
        volScalarField heelMix_;

        //- Pointer to the energy field associated with the vibro-electronic mode [J/kg]
        volScalarField *hevelMix_;

        //- Degrees of freedom associated with the rotational energy mode
        volScalarField zetar_;

        //- Degrees of freedom associated with the vibrational energy mode
        volScalarField zetav_;

        //- Degrees of freedom associated with the electronic energy mode
        volScalarField zetael_;

        //- Switch to downgrade to single vibrational temperature model
        bool downgradeSingleTv_;

        //- Switch to downgrade to single-temperature model
        bool downgradeSingleTemperature_;

        //- Switch to downgrade to single-vibrational mode
        bool downgradeSingleVibMode_;

        //- Run a light version of hy2Foam
        bool hyLight_;


    // Protected Member Functions

        //- Return the trans-rotational enthalpy/internal energy field boundary
        //  types by interrogating the temperature field boundary types
        wordList het2BoundaryTypes();

        //- Return the trans-rotational enthalpy/internal energy field boundary
        //  base types by interrogating the temperature field boundary types
        wordList het2BoundaryBaseTypes();

        //- Return the mixture vibro-electronic enthalpy/internal energy field boundary
        //  types by interrogating the temperature field boundary types
        wordList hevelMix2BoundaryTypes();

        //- Return the mixture vibro-electronic enthalpy/internal energy field boundary
        //  base types by interrogating the temperature field boundary types
        wordList hevelMix2BoundaryBaseTypes();

public:

    //- Runtime type information
    TypeName("multi2Thermo");


    //- Declare run-time constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        multi2Thermo,
        fvMesh,
        (const fvMesh& mesh, const word& phaseName),
        (mesh, phaseName)
    );

    // Constructors

        //- Construct from mesh and phase name
        multi2Thermo
        (
            const fvMesh&,
            const word& phaseName
        );

        //- Construct from mesh and phase name
        multi2Thermo
        (
            const fvMesh&,
            const dictionary&,
            const word& phaseName
        );


    //- Selector
    static autoPtr<multi2Thermo> New
    (
        const fvMesh&,
        const word& phaseName=word::null
    );


    //- Destructor
    virtual ~multi2Thermo();


    // Member functions

        IOdictionary& twoTemperatureDictionary()
        {
            return twoTemperatureDictionary_;
        }

        static const multi2Thermo& lookup2Thermo(const fvPatchScalarField& pf);

        // Downgrade to single vibrational temperature model?
        bool downgradeSingleTv() const;

        // Downgrade to single-temperature model?
        bool downgradeSingleTemperature() const;

        // Downgrade to single-vibrational mode?
        bool downgradeSingleVibMode() const;

        // Run a lighter version of hy2Foam?
        bool hyLight() const;

        // Access to thermodynamic state variables

          //- het
            //- Trans-rotational Enthalpy/Internal energy [J/kg]
            //  Non-const access allowed for transport equations
            virtual volScalarField& het() = 0;

            //- Trans-rotational Enthalpy/Internal energy [J/kg]
            virtual const volScalarField& het() const = 0;

            //- Trans-rotational Enthalpy/Internal energy
            //  for given pressure and temperature [J/kg]
            virtual tmp<volScalarField> het
            (
                const volScalarField& p,
                const volScalarField& Tt
            ) const = 0;

            //- Trans-rotational Enthalpy/Internal energy
            //  for cell-set [J/kg]
            virtual tmp<scalarField> het
            (
                const scalarField& p,
                const scalarField& Tt,
                const labelList& cells
            ) const = 0;

            //- Translational-rotational Enthalpy/Internal energy
            //  for patch [J/kg]
            virtual tmp<scalarField> het
            (
                const scalarField& p,
                const scalarField& Tt,
                const label patchi
            ) const = 0;

          //- hev
            //- Vibrational Enthalpy/Internal energy [J/kg]
            //  Non-const access allowed for transport equations
            virtual volScalarField& hev() = 0;

            //- Vibrational Enthalpy/Internal energy [J/kg]
            virtual const volScalarField& hev() const = 0;

            //- Vibrational Enthalpy/Internal energy
            //  for given pressure and temperature [J/kg]
            virtual tmp<volScalarField> hev
            (
                const volScalarField& p,
                const volScalarField& Tv
            ) const = 0;

            //- Vibrational Enthalpy/Internal energy for cell-set [J/kg]
            virtual tmp<scalarField> hev
            (
                const scalarField& p,
                const scalarField& Tv,
                const labelList& cells
            ) const = 0;

            //- Vibrational Enthalpy/Internal energy for patch [J/kg]
            virtual tmp<scalarField> hev
            (
                const scalarField& p,
                const scalarField& Tv,
                const label patchi
            ) const = 0;

          //- heel
            //- Electronic Enthalpy/Internal energy [J/kg]
            //  Non-const access allowed for transport equations
            virtual volScalarField& heel() = 0;

            //- Electronic Enthalpy/Internal energy [J/kg]
            virtual const volScalarField& heel() const = 0;

            //- Electronic Enthalpy/Internal energy
            //  for given pressure and temperature [J/kg]
            virtual tmp<volScalarField> heel
            (
                const volScalarField& p,
                const volScalarField& Tv
            ) const = 0;

            //- Electronic Enthalpy/Internal energy for cell-set [J/kg]
            virtual tmp<scalarField> heel
            (
                const scalarField& p,
                const scalarField& Tv,
                const labelList& cells
            ) const = 0;

            //- Electronic Enthalpy/Internal energy for patch [J/kg]
            virtual tmp<scalarField> heel
            (
                const scalarField& p,
                const scalarField& Tv,
                const label patchi
            ) const = 0;

          //- hevel
            //- Vib-electronic Enthalpy/Internal energy [J/kg]
            //  Non-const access allowed for transport equations
            virtual volScalarField& hevel() = 0;

            //- Vib-electronic Enthalpy/Internal energy [J/kg]
            virtual const volScalarField& hevel() const = 0;

            //- Vib-electronic Enthalpy/Internal energy
            //  for given pressure and temperature [J/kg]
            virtual tmp<volScalarField> hevel
            (
                const volScalarField& p,
                const volScalarField& Tv
            ) const = 0;

            //- Vib-electronic Enthalpy/Internal energy for cell-set [J/kg]
            virtual tmp<scalarField> hevel
            (
                const scalarField& p,
                const scalarField& Tv,
                const labelList& cells
            ) const = 0;

            //- Vib-electronic Enthalpy/Internal energy for patch [J/kg]
            virtual tmp<scalarField> hevel
            (
                const scalarField& p,
                const scalarField& Tv,
                const label patchi
            ) const = 0;

          /*//- zetar
            //- Degrees of freedom associated with the rotational energy mode
            //  for given pressure and temperature [-]
            virtual tmp<volScalarField> zetar
            (
                const volScalarField& p,
                const volScalarField& Tt,
                const volScalarField& Tv
            ) const = 0;

            //- Degrees of freedom associated with the rotational energy mode
            //  for cell-set [-]
            virtual tmp<scalarField> zetar
            (
                const scalarField& p,
                const scalarField& Tt,
                const scalarField& Tv,
                const labelList& cells
            ) const = 0;

            //- Degrees of freedom assciated with the rotational energy mode
            // for patch [-]
            virtual tmp<scalarField> zetar
            (
                const scalarField& p,
                const scalarField& Tt,
                const scalarField& Tv,
                const label patchi
            ) const = 0;

          //- zetav
            //- Degrees of freedom associated with the vibrational energy mode
            //  for given pressure and temperature [-]
            virtual tmp<volScalarField> zetav
            (
                const volScalarField& p,
                const volScalarField& Tv
            ) const = 0;

            //- Degrees of freedom associated with the vibrational energy mode
            //  for cell-set [-]
            virtual tmp<scalarField> zetav
            (
                const scalarField& p,
                const scalarField& Tv,
                const labelList& cells
            ) const = 0;

            //- Degrees of freedom associated with the vibrational energy mode
            // for patch [-]
            virtual tmp<scalarField> zetav
            (
                const scalarField& p,
                const scalarField& Tv,
                const label patchi
            ) const = 0;
          //- zetael
            //- Degrees of freedom associated with the electronic energy mode
            //  for given pressure and temperature [-]
            virtual tmp<volScalarField> zetael
            (
                const volScalarField& p,
                const volScalarField& Tv
            ) const = 0;

            //- Degrees of freedom associated with the electronic energy mode
            //  for cell-set [-]
            virtual tmp<scalarField> zetael
            (
                const scalarField& p,
                const scalarField& Tv,
                const labelList& cells
            ) const = 0;

            //- Degrees of freedom associated with the electronic energy mode
            // for patch [-]
            virtual tmp<scalarField> zetael
            (
                const scalarField& p,
                const scalarField& Tv,
                const label patchi
            ) const = 0;*/

            //- Heat capacity at constant volume [J/kg/K]
            //  trans-rotational energy mode
            virtual tmp<volScalarField> Cv_t() const = 0;

            //- Heat capacity at constant volume for patch[J/kg/K]
            //  trans-rotational energy mode
            virtual tmp<scalarField> Cv_t
            (
                const scalarField& p,
                const scalarField& Tt,
                const label patchi
            ) const = 0;

            //- Heat capacity at constant volume [J/kg/K]
            //  vib-electronical energy mode
            virtual tmp<volScalarField> Cv_v() const = 0;

            //- Heat capacity at constant volume for patch[J/kg/K]
            //  vib-electronical energy mode
            virtual tmp<scalarField> Cv_v
            (
                const scalarField& p,
                const scalarField& Tv,
                const label patchi
            ) const = 0;

            //- Heat capacity at constant pressure [J/kg/K]
            //  trans-rotational energy mode
            virtual tmp<volScalarField> Cp_t() const = 0;

            //- Heat capacity at constant pressure for patch[J/kg/K]
            //  trans-rotational energy mode
            virtual tmp<scalarField> Cp_t
            (
                const scalarField& p,
                const scalarField& Tt,
                const label patchi
            ) const = 0;

            /*//- Trans-rotational temperature from enthalpy/internal
            //  energy for cell-set
            virtual tmp<scalarField> TtHEt
            (
                const scalarField& h,
                const scalarField& p,
                const scalarField& T0,      // starting temperature
                const labelList& cells
            ) const = 0;

            //- Trans-rotational temperature from enthalpy/internal
            //  energy for patch
            virtual tmp<scalarField> TtHEt
            (
                const scalarField& h,
                const scalarField& p,
                const scalarField& T0,      // starting temperature
                const label patchi
            ) const = 0;

            //- Vibrational temperature from enthalpy/internal
            //  energy for cell-set
            virtual tmp<scalarField> TvHEv
            (
                const scalarField& h,
                const scalarField& p,
                const scalarField& T0,      // starting temperature
                const labelList& cells
            ) const = 0;

            //- Vibrational temperature from enthalpy/internal energy for patch
            virtual tmp<scalarField> TvHEv
            (
                const scalarField& h,
                const scalarField& p,
                const scalarField& T0,      // starting temperature
                const label patchi
            ) const = 0;

            //- Vib-electronic temperature from enthalpy/internal
            //  energy for cell-set
            virtual tmp<scalarField> TvelHEvel
            (
                const scalarField& h,
                const scalarField& p,
                const scalarField& T0,      // starting temperature
                const labelList& cells
            ) const = 0;

            //- Vib-electronic temperature from enthalpy/internal energy for patch
            virtual tmp<scalarField> TvelHEvel
            (
                const scalarField& h,
                const scalarField& p,
                const scalarField& T0,      // starting temperature
                const label patchi
            ) const = 0;*/


        // Fields derived from thermodynamic state variables

            //- Overall temperature of mixture [K]
            virtual const volScalarField& Tov() const;

            //- Overall temperature of mixture for patch [K]
            virtual const scalarField& Tov(const label patchi) const;

            //- Vibrational temperature of mixture [K]
            virtual const volScalarField& Tv() const;

            //- Vibrational temperature of mixture for patch [K]
            virtual const scalarField& Tv(const label patchi) const;

            //- Rotational degrees of freedom of mixture
            virtual const volScalarField& zetar() const;


        // Fields derived from transport state variables

            //- Thermal conductivity of mixture [W/mK], trans-rotational mode
            virtual const volScalarField& kappatr() const = 0;

            //- Thermal conductivity of mixture for patch [W/mK], trans-rotational mode
            virtual const scalarField& kappatr(const label patchi) const = 0;

            //- Return non-const access to the thermal conductivity of mixture [W/mK]
            //  trans-rotational energy mode
            virtual volScalarField& kappatr() = 0;

          //- kappave
            //- Thermal conductivity of mixture [W/mK], vib-electronic mode
            virtual const volScalarField& kappave() const = 0;

            //- Thermal conductivity of mixture for patch [W/mK], vib-electronic mode
            virtual const scalarField& kappave(const label patchi) const = 0;

            //- Return non-const access to the thermal conductivity of mixture [W/mK]
            //  vib-electronic energy mode
            virtual volScalarField& kappave() = 0;

          //- alphatr
            //- Thermal diffusivity of mixture [kg/m/s], trans-rotational mode
            virtual const volScalarField& alphatr() const = 0;

            //- Thermal diffusivity of mixture for patch [kg/m/s], trans-rotational mode
            virtual const scalarField& alphatr(const label patchi) const = 0;

            //- Return non-const access to the thermal diffusivity of mixture [kg/m/s]
            //  trans-rotational energy mode
            virtual volScalarField& alphatr() = 0;

          //- alphave
            //- Thermal diffusivity of mixture [kg/m/s], vib-electronic mode
            virtual const volScalarField& alphave() const = 0;

            //- Thermal diffusivity of mixture for patch [kg/m/s], vib-electronic mode
            virtual const scalarField& alphave(const label patchi) const = 0;

            //- Return non-const access to the thermal diffusivity of mixture [kg/m/s]
            //  vib-electronic energy mode
            virtual volScalarField& alphave() = 0;

            //- Thermal diffusivity of the mixture [kg/m/s]
            //  trans-rotational energy mode
            virtual tmp<volScalarField> alpha_t() const = 0;

            //- Thermal diffusivity of the mixture for patch [kg/m/s]
            //  trans-rotational energy mode
            virtual tmp<scalarField> alpha_t
            (
                const scalarField& p,
                const scalarField& Tt,
                const label patchi
            ) const = 0;

            //- Thermal diffusivity of the mixture [kg/m/s]
            //  vib-electronical energy mode
            virtual tmp<volScalarField> alpha_v() const = 0;

            //- Thermal diffusivity of the mixture for patch [kg/m/s]
            //  vib-electronical energy mode
            virtual tmp<scalarField> alpha_v
            (
                const scalarField& p,
                const scalarField& Tv,
                const label patchi
            ) const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
