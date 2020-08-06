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
    Foam::rho2Thermo

Description
    Basic thermodynamic properties based on density

SourceFiles
    rho2Thermo.C

\*---------------------------------------------------------------------------*/

#ifndef rho2Thermo_H
#define rho2Thermo_H

#include "multi2Thermo.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class rho2Thermo Declaration
\*---------------------------------------------------------------------------*/

class rho2Thermo
:
    public multi2Thermo
{

protected:

    // Protected data

        //- Density field [kg/m^3]
        //  Named 'rho2Thermo' to avoid (potential) conflict with solver density
        volScalarField rho_;

        //- Compressibility [s^2/m^2]
        volScalarField psi_;

        //- Dynamic viscosity [kg/m/s]
        volScalarField mu_;

        //- Thermal conductivity of the mixture [W/mK]
        //  trans-rotational energy mode
        volScalarField kappatr_;

        //- Thermal conductivity of the mixture [W/mK]
        //  vibro-electronic energy mode
        volScalarField kappave_;

        //- Thermal diffusivity of the mixture [kg/ms]
        //  trans-rotational energy mode
        volScalarField alphatr_;

        //- Thermal diffusivity of the mixture [kg/ms]
        //  vibro-electronic energy mode
        volScalarField alphave_;


    // Protected Member Functions

        //- Return the density field boundary types
        wordList rhoBoundaryTypes(); // NEW VINCENT 13/04/2016

        //- Construct as copy (not implemented)
        rho2Thermo(const rho2Thermo&);


public:

    //- Runtime type information
    TypeName("rho2Thermo");


    //- Declare run-time constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        rho2Thermo,
        fvMesh,
        (const fvMesh& mesh, const word& phaseName),
        (mesh, phaseName)
    );


    // Constructors

        //- Construct from mesh and phase name
        rho2Thermo
        (
            const fvMesh&,
            const word& phaseName
        );

        //- Construct from mesh, dictionary and phase name
        rho2Thermo
        (
            const fvMesh&,
            const dictionary&,
            const word& phaseName
        );


    //- Selector
    static autoPtr<rho2Thermo> New
    (
        const fvMesh&,
        const word& phaseName=word::null
    );


    //- Destructor
    virtual ~rho2Thermo();


    // Member functions

        // Fields derived from thermodynamic state variables

          //- rho
            //- Density [kg/m^3]
            virtual tmp<volScalarField> rho() const;

            //- Density for patch [kg/m^3]
            virtual tmp<scalarField> rho(const label patchi) const;

            //- Return non-const access to the local density field [kg/m^3]
            virtual volScalarField& rho();

          //- psi
            //- Compressibility [s^2/m^2]
            virtual const volScalarField& psi() const;

            //- Return non-const access to the local compressibility
            //  field [s^2/m^2]
            virtual volScalarField& getPsi();


        // Access to transport state variables

          //- mu
            //- Dynamic viscosity of mixture [kg/m/s]
            //virtual const volScalarField& mu() const;

            //- Dynamic viscosity of mixture for patch [kg/m/s]
            //virtual const scalarField& mu(const label patchi) const;

            //- Return the dynamic laminar viscosity
            virtual tmp<volScalarField> mu() const;

            //- Return the dynamic laminar viscosity for patch
            virtual tmp<scalarField> mu(const label patchi) const;

            //- Return non-const access to the dynamic viscosity of mixture [kg/m/s]
            virtual volScalarField& mu();

          //- kappatr
            //- Thermal conductivity of mixture [W/mK], trans-rotational mode
            virtual const volScalarField& kappatr() const;

            //- Thermal conductivity of mixture for patch [W/mK], trans-rotational mode
            virtual const scalarField& kappatr(const label patchi) const;

            //- Return non-const access to the thermal conductivity of mixture [W/mK]
            //  trans-rotational energy mode
            virtual volScalarField& kappatr();

          //- kappave
            //- Thermal conductivity of mixture [W/mK], vibro-electronic mode
            virtual const volScalarField& kappave() const;

            //- Thermal conductivity of mixture for patch [W/mK], vibro-electronic mode
            virtual const scalarField& kappave(const label patchi) const;

            //- Return non-const access to the thermal conductivity of mixture [W/mK]
            //  vibro-electronic energy mode
            virtual volScalarField& kappave();

          //- alphatr
            //- Thermal diffusivity of mixture [kg/m/s], trans-rotational mode
            virtual const volScalarField& alphatr() const;

            //- Thermal diffusivity of mixture for patch [kg/m/s], trans-rotational mode
            virtual const scalarField& alphatr(const label patchi) const;

            //- Return non-const access to the thermal diffusivity of mixture [kg/m/s]
            //  trans-rotational energy mode
            virtual volScalarField& alphatr();

          //- alphave
            //- Thermal diffusivity of mixture [kg/m/s], vibro-electronic mode
            virtual const volScalarField& alphave() const;

            //- Thermal diffusivity of mixture for patch [kg/m/s], vibro-electronic mode
            virtual const scalarField& alphave(const label patchi) const;

            //- Return non-const access to the thermal diffusivity of mixture [kg/m/s]
            //  vibro-electronic energy mode
            virtual volScalarField& alphave();

          //- kappaEff
            //- Effective thermal diffusivity for temperature // NEW VINCENT 11/03/2016
            //  of mixture [J/m/s/K]
            virtual tmp<volScalarField> kappaEff(const volScalarField&) const;

            //- Effective thermal diffusivity for temperature
            //  of mixture for patch [J/m/s/K]
            virtual tmp<scalarField> kappaEff
            (
                const scalarField& alphat,
                const label patchi
            ) const;

          //- alphaEff
            //- Effective thermal diffusivity of mixture [kg/m/s]
            virtual tmp<volScalarField> alphaEff
            (
                const volScalarField& alphat
            ) const;

            //- Effective thermal diffusivity of mixture for patch [kg/m/s]
            virtual tmp<scalarField> alphaEff
            (
                const scalarField& alphat,
                const label patchi
            ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
