/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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
    Foam::he2Thermo

Description
    Enthalpy/Internal energy for a mixture using a multi-temperature model

SourceFiles
    he2Thermo.C

\*---------------------------------------------------------------------------*/

#ifndef he2Thermo_H
#define he2Thermo_H

#include "basicMixture.H"
#include "multiThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class he2Thermo Declaration
\*---------------------------------------------------------------------------*/

template<class BasicThermo, class MixtureType>
class he2Thermo
:
    public BasicThermo,
    public MixtureType
{
protected:

    // Protected Member Functions

        // Enthalpy/Internal energy

            //- Correct the enthalpy/internal energy field boundaries
            void hetBoundaryCorrection(volScalarField& h);


private:

    // Private Member Functions

        //- Construct as copy (not implemented)
        he2Thermo(const he2Thermo<BasicThermo, MixtureType>&);


        //- Initialize he2Thermo
        void init();


public:

    // Constructors

        //- Construct from mesh
        he2Thermo
        (
            const fvMesh&,
            const word& phaseName
        );

        //- Construct from mesh and dictionary
        he2Thermo
        (
            const fvMesh&,
            const dictionary&,
            const word& phaseName
        );


    //- Destructor
    virtual ~he2Thermo();


    // Member functions

        //- Return the compostion of the mixture
        virtual typename MixtureType::basicMixtureType&
        composition()
        {
            return *this;
        }

        //- Return the compostion of the mixture
        virtual const typename MixtureType::basicMixtureType&
        composition() const
        {
            return *this;
        }

        // Access to thermodynamic state variables

            // NEW VINCENT ****************************************************
            //- Trans-rotational Enthalpy/Internal energy [J/kg]
            //  Non-const access allowed for transport equations
            virtual volScalarField& het()
            {
                return this->het_;
            }

            //- Trans-rotational Enthalpy/Internal energy [J/kg]
            virtual const volScalarField& het() const
            {
                return this->het_;
            }
            
            //- Vibrational Enthalpy/Internal energy [J/kg]
            //  Non-const access allowed for transport equations
            virtual volScalarField& hev()
            {
                return this->hevMix_;
            }

            //- Vibrational Enthalpy/Internal energy [J/kg]
            virtual const volScalarField& hev() const
            {
                return this->hevMix_;
            }
            
            //- Electronic Enthalpy/Internal energy [J/kg]
            //  Non-const access allowed for transport equations
            virtual volScalarField& heel()
            {
                return this->heelMix_;
            }

            //- Electronic Enthalpy/Internal energy [J/kg]
            virtual const volScalarField& heel() const
            {
                return this->heelMix_;
            }
            
            //- Vib-electronic Enthalpy/Internal energy [J/kg]
            //  Non-const access allowed for transport equations
            virtual volScalarField& hevel()
            {
                return *(this->hevelMix_);
            }

            //- Vib-electronic Enthalpy/Internal energy [J/kg]
            virtual const volScalarField& hevel() const
            {
                return *(this->hevelMix_);
            }
            // END NEW VINCENT ************************************************


        // Fields derived from thermodynamic state variables

            // Old templates with 2 arguments considers thermal equilibrium: Ttr=Tve
            // This modification has been done on 11 Feb 2016
            
            //- Equilibrium heat capacity at constant pressure for patch [J/kg/K]
            virtual tmp<scalarField> Cp
            (
                const scalarField& p,
                const scalarField& T,
                const label patchi
            ) const
            {
                return Cp(p, T, T, patchi);
            };
            
            //- Equilibrium heat capacity at constant volume for patch [J/kg/K]
            virtual tmp<scalarField> Cv
            (
                const scalarField& p,
                const scalarField& T,
                const label patchi
            ) const
            {
                return Cv(p, T, T, patchi);
            };
            
            //- Equilibrium gamma = Cp/Cv for patch []
            virtual tmp<scalarField> gamma
            (
                const scalarField& p,
                const scalarField& T,
                const label patchi
            ) const
            {
                return gamma(p, T, T, patchi);
            };
            
            // NEW VINCENT ****************************************************
          //- het
            //- Trans-rotational Enthalpy/Internal energy
            //  for given pressure and temperature [J/kg]
            virtual tmp<volScalarField> het
            (
                const volScalarField& p,
                const volScalarField& Tt
            ) const;

            //- Trans-rotational Enthalpy/Internal 
            //  energy for cell-set [J/kg]
            virtual tmp<scalarField> het
            (
                const scalarField& p,
                const scalarField& Tt,
                const labelList& cells
            ) const;

            //- Trans-rotational Enthalpy/Internal 
            //  energy for patch [J/kg]
            virtual tmp<scalarField> het
            (
                const scalarField& p,
                const scalarField& Tt,
                const label patchi
            ) const;
            
          //- hev
            //- Vibrational Enthalpy/Internal energy
            //  for given pressure and temperature [J/kg]
            virtual tmp<volScalarField> hev
            (
                const volScalarField& p,
                const volScalarField& Tv
            ) const;

            //- Vibrational Enthalpy/Internal energy for cell-set [J/kg]
            virtual tmp<scalarField> hev
            (
                const scalarField& p,
                const scalarField& Tv,
                const labelList& cells
            ) const;

            //- Vibrational Enthalpy/Internal energy for patch [J/kg]
            virtual tmp<scalarField> hev
            (
                const scalarField& p,
                const scalarField& Tv,
                const label patchi
            ) const;
            
          //- heel  
            //- Electronic Enthalpy/Internal energy
            //  for given pressure and temperature [J/kg]
            virtual tmp<volScalarField> heel
            (
                const volScalarField& p,
                const volScalarField& Tv
            ) const;

            //- Electronic Enthalpy/Internal energy for cell-set [J/kg]
            virtual tmp<scalarField> heel
            (
                const scalarField& p,
                const scalarField& Tv,
                const labelList& cells
            ) const;
            
            //- Electronic Enthalpy/Internal energy for patch [J/kg]
            virtual tmp<scalarField> heel
            (
                const scalarField& p,
                const scalarField& Tv,
                const label patchi
            ) const;
            
          //- hevel  
            //- Vib-electronic Enthalpy/Internal energy
            //  for given pressure and temperature [J/kg]
            virtual tmp<volScalarField> hevel
            (
                const volScalarField& p,
                const volScalarField& Tv
            ) const;

            //- Vib-electronic Enthalpy/Internal energy for cell-set [J/kg]
            virtual tmp<scalarField> hevel
            (
                const scalarField& p,
                const scalarField& Tv,
                const labelList& cells
            ) const;
            
            //- Vib-electronic Enthalpy/Internal energy for patch [J/kg]
            virtual tmp<scalarField> hevel
            (
                const scalarField& p,
                const scalarField& Tv,
                const label patchi
            ) const;
            
          //- zetar  
            //- Degrees of freedom associated with the rotational energy mode
            //  for given pressure and temperature [-]
            virtual tmp<volScalarField> zetar
            (
                const volScalarField& p,
                const volScalarField& Tt,
                const volScalarField& Tv
            ) const;

            //- Degrees of freedom associated with the rotational energy mode
            //  for cell-set [-]
            virtual tmp<scalarField> zetar
            (
                const scalarField& p,
                const scalarField& Tt,
                const scalarField& Tv,
                const labelList& cells
            ) const;

            //- Degree of freedom assciated with the rotational energy mode
            // for patch [-]
            virtual tmp<scalarField> zetar
            (
                const scalarField& p,
                const scalarField& Tt,
                const scalarField& Tv,
                const label patchi
            ) const; 
            
          //- zetav 
            //- Degrees of freedom associated with the vibrational energy mode
            //  for given pressure and temperature [-]
            virtual tmp<volScalarField> zetav
            (
                const volScalarField& p,
                const volScalarField& Tv
            ) const;

            //- Degrees of freedom associated with the vibrational energy mode
            //  for cell-set [-]
            virtual tmp<scalarField> zetav
            (
                const scalarField& p,
                const scalarField& Tv,
                const labelList& cells
            ) const;

            //- Degrees of freedom associated with the vibrational energy mode
            // for patch [-]
            virtual tmp<scalarField> zetav
            (
                const scalarField& p,
                const scalarField& Tv,
                const label patchi
            ) const;
            
          //- zetael 
            //- Degrees of freedom associated with the electronic energy mode
            //  for given pressure and temperature [-]
            virtual tmp<volScalarField> zetael
            (
                const volScalarField& p,
                const volScalarField& Tv
            ) const;

            //- Degrees of freedom associated with the electronic energy mode
            //  for cell-set [-]
            virtual tmp<scalarField> zetael
            (
                const scalarField& p,
                const scalarField& Tv,
                const labelList& cells
            ) const;

            //- Degrees of freedom associated with the electronic energy mode
            // for patch [-]
            virtual tmp<scalarField> zetael
            (
                const scalarField& p,
                const scalarField& Tv,
                const label patchi
            ) const;               
            // END NEW VINCENT ************************************************

            //- Chemical enthalpy [J/kg]
            virtual tmp<volScalarField> hc() const;

            // NEW VINCENT ****************************************************
            //- Trans-rotational temperature from enthalpy/internal 
            //  energy for cell-set
            virtual tmp<scalarField> TtHEt
            (
                const scalarField& h,
                const scalarField& p,
                const scalarField& T0,      // starting temperature
                const labelList& cells
            ) const;

            //- Trans-rotational temperature from enthalpy/internal
            //  energy for cell-set
            virtual tmp<scalarField> TtHEt
            (
                const scalarField& h,
                const scalarField& p,
                const scalarField& T0,      // starting temperature
                const label patchi
            ) const;
            
            //- Vibrational temperature from enthalpy/internal 
            //  energy for cell-set
            virtual tmp<scalarField> TvHEv
            (
                const scalarField& h,
                const scalarField& p,
                const scalarField& T0,      // starting temperature
                const labelList& cells
            ) const;

            //- Vibrational temperature from enthalpy/internal energy for patch
            virtual tmp<scalarField> TvHEv
            (
                const scalarField& h,
                const scalarField& p,
                const scalarField& T0,      // starting temperature
                const label patchi
            ) const;
            
            //- Vib-electronic temperature from enthalpy/internal 
            //  energy for cell-set
            virtual tmp<scalarField> TvelHEvel
            (
                const scalarField& h,
                const scalarField& p,
                const scalarField& T0,      // starting temperature
                const labelList& cells
            ) const;

            //- Vib-electronic temperature from enthalpy/internal energy for patch
            virtual tmp<scalarField> TvelHEvel
            (
                const scalarField& h,
                const scalarField& p,
                const scalarField& T0,      // starting temperature
                const label patchi
            ) const;
            // END NEW VINCENT ************************************************
            
            //- Heat capacity at constant pressure for patch [J/kg/K]
            virtual tmp<scalarField> Cp
            (
                const scalarField& p,
                const scalarField& Tt,
                const scalarField& Tv,
                const label patchi
            ) const;

            //- Heat capacity at constant pressure [J/kg/K]
            virtual tmp<volScalarField> Cp() const;

            //- Heat capacity at constant volume for patch [J/kg/K]
            virtual tmp<scalarField> Cv
            (
                const scalarField& p,
                const scalarField& Tt,
                const scalarField& Tv,
                const label patchi
            ) const;

            //- Heat capacity at constant volume [J/kg/K]
            virtual tmp<volScalarField> Cv() const;
            
            // NEW VINCENT **********************************************
            //- Heat capacity at constant volume [J/kg/K]
            //  trans-rotational energy mode
            virtual tmp<volScalarField> Cv_t() const;

            //- Heat capacity at constant volume for patch[J/kg/K]
            //  trans-rotational energy mode
            virtual tmp<scalarField> Cv_t
            (
                const scalarField& p,
                const scalarField& Tt,
                const label patchi
            ) const;
            
            //- Heat capacity at constant volume [J/kg/K]
            //  vibrational energy mode
            virtual tmp<volScalarField> Cv_v() const;
            
            //- Heat capacity at constant volume for patch[J/kg/K]
            //  vibrational energy mode
            virtual tmp<scalarField> Cv_v
            (
                const scalarField& p,
                const scalarField& Tv,
                const label patchi
            ) const;
            
            //- Heat capacity at constant pressure [J/kg/K]
            //  trans-rotational energy mode
            virtual tmp<volScalarField> Cp_t() const;

            //- Heat capacity at constant pressure for patch[J/kg/K]
            //  trans-rotational energy mode
            virtual tmp<scalarField> Cp_t
            (
                const scalarField& p,
                const scalarField& Tt,
                const label patchi
            ) const;
            // END NEW VINCENT ******************************************

            //- gamma = Cp/Cv []
            virtual tmp<volScalarField> gamma() const;

            //- gamma = Cp/Cv for patch []
            virtual tmp<scalarField> gamma
            (
                const scalarField& p,
                const scalarField& Tt,
                const scalarField& Tv,
                const label patchi
            ) const;


        // Fields derived from transport state variables

            // NEW VINCENT **********************************************
            //- Thermal diffusivity of the mixture [kg/m/s]
            //  trans-rotational energy mode
            virtual tmp<volScalarField> alpha_t() const;

            //- Thermal diffusivity of the mixture for patch [kg/m/s]
            //  trans-rotational energy mode
            virtual tmp<scalarField> alpha_t
            (
                const scalarField& p,
                const scalarField& Tt,
                const label patchi
            ) const;
            
            //- Thermal diffusivity of the mixture [kg/m/s]
            //  vibrational energy mode
            virtual tmp<volScalarField> alpha_v() const;

            //- Thermal diffusivity of the mixture for patch [kg/m/s]
            //  vibrational energy mode
            virtual tmp<scalarField> alpha_v
            (
                const scalarField& p,
                const scalarField& Tv,
                const label patchi
            ) const;
             // END BRAND NEW VINCENT ******************************************

        //- Read thermophysical properties dictionary
        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

#ifdef NoRepository
#   include "he2Thermo.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
