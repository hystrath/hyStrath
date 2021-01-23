/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016-2021 hyStrath
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
    Foam::Specie2Mixture

Description
    Foam::Specie2Mixture

SourceFiles
    Specie2Mixture.C

\*---------------------------------------------------------------------------*/

#ifndef Specie2Mixture_H
#define Specie2Mixture_H

#include "scalar.H"

#include "volFields.H"
#include "DynamicList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class dictionary;
class fvMesh;


/*---------------------------------------------------------------------------*\
                     Class Specie2Mixture Declaration
\*---------------------------------------------------------------------------*/

template<class MixtureType>
class Specie2Mixture
:
    public MixtureType
{

public:

    // Constructors

        //- Construct from dictionary and mesh
        Specie2Mixture(const dictionary&, const fvMesh&);


    //- Destructor
    virtual ~Specie2Mixture()
    {}


    // Member functions

        // Per specie properties

            //- Number of moles []
            virtual scalar nMoles(const label speciei) const;

            //- Molecular weight [kg/kmol]
            virtual scalar W(const label speciei) const;

            //- Molecular weight [kg/kmol]
            virtual scalar W(const word& specieName) const;
            
            //- Species gas constant [J/kg/K]
            virtual scalar R(const label speciei) const;

            //- Particle type
            virtual label particleType(const label speciei) const;

            //- Particle charge
            virtual label particleCharge(const label speciei) const;

            //- Diameter of the particle
            virtual scalar diameter(const label speciei) const;

            //- Vibrational list with degeneracy and characteristic vibrational
            //  temperature [K] for each vib. mode
            virtual DynamicList<scalar> vibrationalList
            (
                const label speciei
            ) const;

            //- Dissociation potential of the particle [J/kg]
            virtual scalar dissociationPotential(const label speciei) const;

            //- Temperature-viscosity coefficient
            virtual scalar omega(const label speciei) const;

            //- Number of characteristic vibrational temperatures
            virtual label noVibrationalTemp(const label speciei) const;

            //- Number of electronic energy levels
            virtual label noElectronicLevels(const label speciei) const;

            //- Fraction of the first ionisation energy of a given species [J/kg]
            virtual scalar iHat(const label speciei) const;


        // Per specie thermo properties

            // Backward compatibility
            //- Heat capacity at constant pressure [J/(kg K)]
            virtual scalar Cp
            (
                const label speciei,
                const scalar p,
                const scalar T
            ) const;

            //- Heat capacity at constant volume [J/(kg K)]
            virtual scalar Cv
            (
                const label speciei,
                const scalar p,
                const scalar T
            ) const;

            //- Absolute enthalpy [J/kg]
            virtual scalar Ha
            (
                const label speciei,
                const scalar p,
                const scalar T
            ) const;

            //- Sensible enthalpy [J/kg]
            virtual scalar Hs
            (
                const label speciei,
                const scalar p,
                const scalar T
            ) const;
            // End backward compatibility

            //- Chemical enthalpy [J/kg]
            virtual scalar Hc(const label speciei) const;

            // Backward compatibility
            //- Entropy [J/(kg K)]
            virtual scalar S
            (
                const label speciei,
                const scalar p,
                const scalar T
            ) const;

            //- Sensible internal energy [J/kg]
            virtual scalar Es
            (
                const label speciei,
                const scalar p,
                const scalar T
            ) const;

            //- Gibbs free energy [J/kg]
            virtual scalar G
            (
                const label speciei,
                const scalar p,
                const scalar T
            ) const;

            //- Helmholtz free energy [J/kg]
            virtual scalar A
            (
                const label speciei,
                const scalar p,
                const scalar T
            ) const;
            // End backward compatibility

            //- Heat capacity at constant pressure [J/(kg K)]
            virtual scalar Cp
            (
                const label speciei,
                const scalar p,
                const scalar T,
                const scalar Tve
            ) const;

            //- Heat capacity at constant volume [J/(kg K)]
            virtual scalar Cv
            (
                const label speciei,
                const scalar p,
                const scalar T,
                const scalar Tve
            ) const;

            //- Heat capacity at constant volume [J/(kg K)]
            //  trans-rotational energy mode
            virtual scalar Cv_t
            (
                const label speciei,
                const scalar p,
                const scalar T
            ) const;

            //- Heat capacity at constant volume [J/(kg K)]
            //  vibrational energy mode
            virtual scalar Cv_v
            (
                const label speciei,
                const scalar p,
                const scalar Tve
            ) const;

            //- Heat capacity at constant volume [J/(kg K)]
            //  electronic energy mode
            virtual scalar Cv_el
            (
                const label speciei,
                const scalar p,
                const scalar Tve
            ) const;

            //- Heat capacity at constant pressure [J/(kg K)]
            //  Vibro-electronic energy mode
            virtual scalar Cv_vel
            (
                const label speciei,
                const scalar p,
                const scalar Tve
            ) const;

            //- Heat capacity at constant pressure [J/(kg K)]
            //  trans-rotational energy mode
            virtual scalar Cp_t
            (
                const label speciei,
                const scalar p,
                const scalar T
            ) const;

            //- Heat capacity at constant pressure [J/(kg K)]
            //  vibrational energy mode
            virtual scalar Cp_v
            (
                const label speciei,
                const scalar p,
                const scalar Tve
            ) const;

            //- Heat capacity at constant pressure [J/(kg K)]
            //  electronic energy mode
            virtual scalar Cp_el
            (
                const label speciei,
                const scalar p,
                const scalar Tve
            ) const;

            //- Heat capacity at constant pressure [J/(kg K)]
            //  Vibro-electronic energy mode
            virtual scalar Cp_vel
            (
                const label speciei,
                const scalar p,
                const scalar Tve
            ) const;

            //- Absolute enthalpy [J/kg]
            virtual scalar Ha
            (
                const label speciei,
                const scalar p,
                const scalar T,
                const scalar Tve
            ) const;

            //- Sensible enthalpy [J/kg]
            virtual scalar Hs
            (
                const label speciei,
                const scalar p,
                const scalar T,
                const scalar Tve
            ) const;

            //- Sensible enthalpy, trans-rotational [J/kg]
            virtual scalar Hts
            (
                const label speciei,
                const scalar p,
                const scalar T
            ) const;

            //- Sensible enthalpy, vibrational [J/kg]
            virtual scalar Hvs
            (
                const label speciei,
                const scalar p,
                const scalar Tve
            ) const;

            //- Sensible enthalpy, electronic [J/kg]
            virtual scalar Hels
            (
                const label speciei,
                const scalar p,
                const scalar Tve
            ) const;

            //- Sensible enthalpy, Vibro-electronic [J/kg]
            virtual scalar Hvels
            (
                const label speciei,
                const scalar p,
                const scalar Tve
            ) const;

            //- Entropy [J/(kg K)]
            virtual scalar S
            (
                const label speciei,
                const scalar p,
                const scalar T,
                const scalar Tve
            ) const;

            //- Sensible internal energy [J/kg]
            virtual scalar Es
            (
                const label speciei,
                const scalar p,
                const scalar T,
                const scalar Tve
            ) const;

            //- Sensible internal energy, trans-rot. [J/kg]
            virtual scalar Ets
            (
                const label speciei,
                const scalar p,
                const scalar T
            ) const;

            //- Sensible internal energy, vibration. [J/kg]
            virtual scalar Evs
            (
                const label speciei,
                const scalar p,
                const scalar Tve
            ) const;

            //- Sensible internal energy, electronic [J/kg]
            virtual scalar Eels
            (
                const label speciei,
                const scalar p,
                const scalar Tve
            ) const;

            //- Sensible internal energy, Vibro-electronic [J/kg]
            virtual scalar Evels
            (
                const label speciei,
                const scalar p,
                const scalar Tve
            ) const;

             //- Trans-rotational Enthalpy/Internal energy [J/kg]
            virtual scalar HEt
            (
                const label speciei,
                const scalar p,
                const scalar T
            ) const;

            //- Vibrational Enthalpy/Internal energy [J/kg]
            virtual scalar HEv
            (
                const label speciei,
                const scalar p,
                const scalar Tve
            ) const;

            //- Electronic Enthalpy/Internal energy [J/kg]
            virtual scalar HEel
            (
                const label speciei,
                const scalar p,
                const scalar Tve
            ) const;

            //- Vibro-electronic Enthalpy/Internal energy [J/kg]
            virtual scalar HEvel
            (
                const label speciei,
                const scalar p,
                const scalar Tve
            ) const;

            //- Vibro-electronic Enthalpy/Internal energy [J/kg]
            //  for the vibrational mode m TODO ABORTIVE WORK
            /*virtual scalar HEvel_mode
            (
                const label speciei,
                const label mode,
                const scalar p,
                const scalar Tve
            ) const;*/

            //- Heat capacity at constant volume for patch[J/kg/K]
            //  Trans-rotational energy mode
            virtual tmp<scalarField> Cv_t
            (
                const label speciei,
                const scalarField& p,
                const scalarField& T,
                const label patchi
            ) const;
            
            //- Heat capacity at constant volume for patch[J/kg/K]
            //  Vibro-electronic energy mode
            virtual tmp<scalarField> Cv_vel
            (
                const label speciei,
                const scalarField& p,
                const scalarField& Tve,
                const label patchi
            ) const;

            //- Trans-rotational Enthalpy/Internal energy for cell-set [J/kg]
            virtual tmp<scalarField> het
            (
                const label speciei,
                const scalarField& p,
                const scalarField& T,
                const labelList& cells
            ) const;

            //- Trans-rotational Enthalpy/Internal energy for patch [J/kg]
            virtual tmp<scalarField> het
            (
                const label speciei,
                const scalarField& p,
                const scalarField& T,
                const label patchi
            ) const;
            
            //- Vibro-electronic Enthalpy/Internal energy for cell-set [J/kg]
            virtual tmp<scalarField> hevel
            (
                const label speciei,
                const scalarField& p,
                const scalarField& Tve,
                const labelList& cells
            ) const;

            //- Vibro-electronic Enthalpy/Internal energy for patch [J/kg]
            virtual tmp<scalarField> hevel
            (
                const label speciei,
                const scalarField& p,
                const scalarField& Tve,
                const label patchi
            ) const;

            //- Trans-rotational temperature from sensible internal energy
            //  given an initial T0
            virtual scalar TtHEt
            (
                const label speciei,
                const scalar het,
                const scalar p,
                const scalar T
            ) const;

            //- Degrees of freedom associated with the rot mode
            virtual scalar zetar
            (
                const label speciei,
                const scalar p,
                const scalar T,
                const scalar Tve
            ) const;

            //- Degrees of freedom associated with the vib mode
            virtual scalar zetav
            (
                const label speciei,
                const scalar p,
                const scalar Tve
            ) const;

            //- Degrees of freedom associated with the vib mode
            //  for the vibrational mode m TODO ABORTIVE WORK
            /*virtual scalar zetav_mode
            (
                const label speciei,
                const label mode,
                const scalar p,
                const scalar Tve
            ) const;*/

            //- Degrees of freedom associated with the electronic energy mode
            virtual scalar zetael
            (
                const label speciei,
                const scalar p,
                const scalar Tve
            ) const;

            //- Gibbs free energy [J/kg]
            virtual scalar G
            (
                const label speciei,
                const scalar p,
                const scalar T,
                const scalar Tve
            ) const;

            //- Helmholtz free energy [J/kg]
            virtual scalar A
            (
                const label speciei,
                const scalar p,
                const scalar T,
                const scalar Tve
            ) const;


        // Per specie transport properties

            //- Compressibility rho/p [s^2/m^2]
            virtual scalar psi
            (
                const label speciei,
                const scalar p,
                const scalar T
            ) const;

            //- Dynamic viscosity [kg/m/s]
            virtual scalar mu
            (
                const label speciei,
                const scalar p,
                const scalar T
            ) const;

            //- Thermal conductivity [W/m/K]
            //  trans-rotational energy mode
            virtual scalar kappatr
            (
                const label speciei,
                const scalar p,
                const scalar T
            ) const;

            //- Thermal conductivity [W/m/K]
            //  vibrational-electron-electronic energy mode
            virtual scalar kappave
            (
                const label speciei,
                const scalar p,
                const scalar T,
                const scalar Tve = 0
            ) const;
            
            //- Thermal conductivity [W/m/K]
            virtual scalar kappa
            (
                const label speciei,
                const scalar p,
                const scalar T,
                const scalar Tve = 0
            ) const;

            //- Thermal diffusivity [kg/m/s]
            //  trans-rotational energy mode
            virtual scalar alphatr
            (
                const label speciei,
                const scalar p,
                const scalar T
            ) const;

            //- Thermal diffusivity [kg/m/s]
            //  Vibro-electronic energy mode
            virtual scalar alphave
            (
                const label speciei,
                const scalar p,
                const scalar T,
                const scalar Tve = 0
            ) const;
            
            //- Thermal diffusivity [kg/m/s]
            virtual scalar alpha
            (
                const label speciei,
                const scalar p,
                const scalar T,
                const scalar Tve = 0
            ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

#ifdef NoRepository
#   include "Specie2Mixture.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
