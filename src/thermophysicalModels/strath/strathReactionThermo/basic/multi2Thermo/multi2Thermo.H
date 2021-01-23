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

#include "fluid2Thermo.H"

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

        //- Pointer to the energy field associated with the trans-rotational
        //  energy mode [J/kg]
        volScalarField *het_;

        //- Energy field associated with the vibrational energy mode [J/kg]
        volScalarField hevMix_;

        //- Energy field associated with the electronic energy mode [J/kg]
        volScalarField heelMix_;

        //- Pointer to the energy field associated with the vibro-electronic
        //  energy mode [J/kg]
        volScalarField *hevelMix_;
        
        //- Mixture heat capacity at constant volume, trans-rotational energy
        //  mode [J/kg/K]
        volScalarField CvtrMix_;
        
        //- Mixture heat capacity at constant pressure, trans-rotational energy
        //  mode [J/kg/K]
        volScalarField CptrMix_;
        
        //- Mixture heat capacity at constant volume, vibro-electronic energy
        //  mode [J/kg/K]
        volScalarField CvvelMix_;
        
        //- Mixture heat capacity at constant pressure, vibro-electronic energy
        //  mode [J/kg/K]
        volScalarField CpvelMix_;
        
        //- Mixture gas constant [J/kg/K]
        volScalarField RMix_;
        
        //- Mixture heat capacity at constant volume [J/kg/K]
        volScalarField CvMix_;
        
        //- Mixture heat capacity at constant pressure [J/kg/K]
        volScalarField CpMix_;
        
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

        //- Return the mixture vibro-electronic enthalpy/internal energy field
        //  boundary types by interrogating the temperature field boundary types
        wordList hevelMix2BoundaryTypes();

        //- Return the mixture vibro-electronic enthalpy/internal energy field
        //  boundary base types by interrogating the temperature field boundary
        //  types
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

            //- Trans-rotational enthalpy/internal energy [J/kg]
            //  Non-const access allowed for transport equations
            virtual volScalarField& het() = 0;

            //- Trans-rotational enthalpy/internal energy [J/kg]
            virtual const volScalarField& het() const = 0;

            //- Trans-rotational enthalpy/internal energy
            //  for given pressure and temperature [J/kg]
            virtual tmp<volScalarField> het
            (
                const volScalarField& p,
                const volScalarField& T
            ) const = 0;

            //- Trans-rotational enthalpy/internal energy for cell-set [J/kg]
            virtual tmp<scalarField> het
            (
                const scalarField& p,
                const scalarField& T,
                const labelList& cells
            ) const = 0;

            //- Trans-rotational enthalpy/internal energy for patch [J/kg]
            virtual tmp<scalarField> het
            (
                const scalarField& p,
                const scalarField& T,
                const label patchi
            ) const = 0;

            //- Vibrational enthalpy/internal energy [J/kg]
            //  Non-const access allowed for transport equations
            virtual volScalarField& hev() = 0;

            //- Vibrational enthalpy/internal energy [J/kg]
            virtual const volScalarField& hev() const = 0;

            //- Electronic enthalpy/internal energy [J/kg]
            //  Non-const access allowed for transport equations
            virtual volScalarField& heel() = 0;

            //- Electronic enthalpy/internal energy [J/kg]
            virtual const volScalarField& heel() const = 0;

            //- Vibro-electronic enthalpy/internal energy [J/kg]
            //  Non-const access allowed for transport equations
            virtual volScalarField& hevel() = 0;

            //- Vibro-electronic enthalpy/internal energy [J/kg]
            virtual const volScalarField& hevel() const = 0;

            //- Heat capacity at constant volume [J/kg/K]
            //  trans-rotational energy mode
            virtual tmp<volScalarField> Cv_t() const = 0;

            //- Heat capacity at constant volume for patch[J/kg/K]
            //  trans-rotational energy mode
            virtual tmp<scalarField> Cv_t
            (
                const scalarField& p,
                const scalarField& T,
                const label patchi
            ) const = 0;

            //- Heat capacity at constant volume [J/kg/K]
            //  vibro-electronic energy mode
            virtual tmp<volScalarField> Cv_v() const = 0;

            //- Heat capacity at constant pressure [J/kg/K]
            //  trans-rotational energy mode
            virtual tmp<volScalarField> Cp_t() const;

            //- Heat capacity at constant pressure for patch [J/kg/K]
            //  trans-rotational energy mode
            virtual tmp<scalarField> Cp_t
            (
                const scalarField& p,
                const scalarField& T,
                const label patch
            ) const = 0;
            
            //- Heat capacity at constant volume of mixture [J/kg/K]
            //  trans-rotational energy mode
            virtual const volScalarField& CvtrMix() const;
            
            //- Non-const access to heat capacity at constant volume of
            //  mixture [J/kg/K], trans-rotational energy mode
            virtual volScalarField& CvtrMix();
            
            //- Heat capacity at constant pressure of mixture [J/kg/K]
            //  trans-rotational energy mode
            virtual const volScalarField& CptrMix() const;
            
            //- Non-const access to heat capacity at constant pressure of
            //  mixture [J/kg/K], trans-rotational energy mode
            virtual volScalarField& CptrMix();
            
            //- Heat capacity at constant volume of mixture [J/kg/K]
            //  vibro-electronic energy mode
            virtual const volScalarField& CvvelMix() const;
            
            //- Non-const access to heat capacity at constant volume of
            //  mixture [J/kg/K], vibro-electronic energy mode
            virtual volScalarField& CvvelMix();
            
            //- Heat capacity at constant pressure of mixture [J/kg/K]
            //  vibro-electronic energy mode
            virtual const volScalarField& CpvelMix() const;
            
            //- Non-const access to heat capacity at constant pressure of
            //  mixture [J/kg/K], vibro-electronic energy mode
            virtual volScalarField& CpvelMix();
            
            //- Heat capacity at constant volume of mixture [J/kg/K]
            virtual const volScalarField& CvMix() const;
            
            //- Non-const access to heat capacity at constant volume of
            //  mixture [J/kg/K]
            virtual volScalarField& CvMix();
            
            //- Heat capacity at constant pressure of mixture [J/kg/K]
            virtual const volScalarField& CpMix() const;
            
            //- Heat capacity at constant pressure of mixture for patch [J/kg/K]
            virtual const scalarField& CpMix(const label patchi) const;
            
            //- Non-const access to heat capacity at constant pressure of
            //  mixture [J/kg/K]
            virtual volScalarField& CpMix();


        // Fields derived from thermodynamic state variables

            //- Overall temperature of mixture [K]
            virtual const volScalarField& Tov() const;

            //- Overall temperature of mixture for patch [K]
            virtual const scalarField& Tov(const label patchi) const;
            
            //- Non-const access to overall temperature of mixture [K]
            virtual volScalarField& Tov();

            //- Vibro-electronic temperature of mixture [K]
            virtual const volScalarField& Tv() const;
            
            //- Vibro-electronic temperature of mixture for patch [K]
            virtual const scalarField& Tv(const label patchi) const;
            
            //- Non-const access to vibro-electronic temperature of mixture [K]
            virtual volScalarField& Tv();
            
            //- Mixture gas constant [J/kg/K]
            virtual const volScalarField& RMix() const;
            
            //- Mixture gas constant for patch [J/kg/K]
            virtual const scalarField& RMix(const label patchi) const;
            
            //- Non-const access to mixture gas constant [J/kg/K]
            virtual volScalarField& RMix();

            //- Rotational degrees of freedom of mixture
            virtual const volScalarField& zetar() const;
            
            //- Non-const access to rotational degrees of freedom of mixture
            virtual volScalarField& zetar();
            
            //- Vibrational degrees of freedom of mixture
            virtual const volScalarField& zetav() const;
            
            //- Non-const access to vibrational degrees of freedom of mixture
            virtual volScalarField& zetav();
            
            //- Electronic degrees of freedom of mixture
            virtual const volScalarField& zetael() const;
            
            //- Non-const access to electronic degrees of freedom of mixture
            virtual volScalarField& zetael();


        // Fields derived from transport state variables

            //- Thermal conductivity of mixture [W/mK]
            //  trans-rotational energy mode
            virtual const volScalarField& kappatr() const = 0;

            //- Thermal conductivity of mixture for patch [W/mK]
            //  trans-rotational energy mode
            virtual const scalarField& kappatr(const label patchi) const = 0;

            //- Return non-const access to the thermal conductivity of
            //  mixture [W/mK], trans-rotational energy mode
            virtual volScalarField& kappatr() = 0;

            //- Thermal conductivity of mixture [W/mK]
            //  vibro-electronic energy mode
            virtual const volScalarField& kappave() const = 0;

            //- Thermal conductivity of mixture for patch [W/mK]
            //  vibro-electronic energy mode
            virtual const scalarField& kappave(const label patchi) const = 0;

            //- Return non-const access to the thermal conductivity of
            //  mixture [W/mK], vibro-electronic energy mode
            virtual volScalarField& kappave() = 0;
            
            //- Thermal conductivity of mixture [W/mK]
            virtual const volScalarField& kappa() const = 0;

            //- Thermal conductivity of mixture for patch [W/mK]
            virtual const scalarField& kappa(const label patchi) const = 0;

            //- Return non-const access to the thermal conductivity of
            virtual volScalarField& kappa() = 0;

            //- Thermal diffusivity of mixture [kg/m/s]
            //  trans-rotational energy mode
            virtual const volScalarField& alphatr() const = 0;

            //- Thermal diffusivity of mixture for patch [kg/m/s]
            //  trans-rotational energy mode
            virtual const scalarField& alphatr(const label patchi) const = 0;

            //- Return non-const access to the thermal diffusivity of
            //  mixture [kg/m/s], trans-rotational energy mode
            virtual volScalarField& alphatr() = 0;

            //- Thermal diffusivity of mixture [kg/m/s]
            //  vibro-electronic energy mode
            virtual const volScalarField& alphave() const = 0;

            //- Thermal diffusivity of mixture for patch [kg/m/s]
            //  vibro-electronic energy mode
            virtual const scalarField& alphave(const label patchi) const = 0;

            //- Return non-const access to the thermal diffusivity of
            //  mixture [kg/m/s], vibro-electronic energy mode
            virtual volScalarField& alphave() = 0;
            
            //- Thermal diffusivity of mixture [kg/m/s]
            virtual const volScalarField& alpha() const = 0;

            //- Thermal diffusivity of mixture for patch [kg/m/s]
            virtual const scalarField& alpha(const label patchi) const = 0;

            //- Return non-const access to the thermal diffusivity of
            //  mixture [kg/m/s]
            virtual volScalarField& alpha() = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
