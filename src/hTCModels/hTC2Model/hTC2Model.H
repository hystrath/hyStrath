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
    Foam::hTC2Model

Description
    Base class for hTC models

SourceFiles
    hTC2Model.C

\*---------------------------------------------------------------------------*/

#ifndef hTC2Model_H
#define hTC2Model_H

#include "IOdictionary.H"
#include "turbulentFluidThermoModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class hTC2Model Declaration
\*---------------------------------------------------------------------------*/

class hTC2Model
:
    public IOdictionary
{

protected:

    // Protected data

        //- Reference to the turbulence model
        compressibleTurbulenceModel* turbulencePtr_;

        //- Reference to the mesh database
        const fvMesh& mesh_;

        //- Active
        Switch active_;

        //- Dictionary of the model
        dictionary coeffs_;

        //- Model name
        const word modelType_;


private:

    // Private Member Functions

        //- Disallow copy construct
        hTC2Model(const hTC2Model&);

        //- Disallow default bitwise assignment
        void operator=(const hTC2Model&);


public:

    //- Runtime type information
    TypeName("hTC2Model");


    // Constructors


        //- Construct from components
        hTC2Model(const word& modelType, const fvMesh& mesh);


    //- Destructor
    virtual ~hTC2Model();


    // Member Functions

        // Access

            //- Return const access to the mesh database
            inline const fvMesh& mesh() const;

            //- Return const access to phi
            inline const surfaceScalarField& phi() const;

            //- Return const access to rho
            virtual tmp<volScalarField> rho() const = 0;

            //- Return access to turbulence
            inline const compressibleTurbulenceModel& turbulence() const;

            //- Set turbulence
            inline void setTurbulence(compressibleTurbulenceModel& turbModel);

            //- Is hTC active?
            inline const Switch& active() const;

            //- Return const dictionary of the model
            inline const dictionary& coeffs() const;


    // Evolution

        //- Correct hTC rate
        virtual void correct() = 0;

        //- Fuel consumption rate matrix, i.e. source term for fuel equation
        virtual tmp<fvScalarMatrix> R(volScalarField& Y) const = 0;

        //- Heat release rate calculated from fuel consumption rate matrix
        virtual tmp<volScalarField> dQ() const = 0;

        //-  Return source for enthalpy equation [kg/m/s3]
        virtual tmp<volScalarField> Sh() const;

        //-  Return chemistry-vibrational coupling source term
        //  in the vibrational equation [kg/m/s3]
        virtual tmp<volScalarField> Scv() const;
        
        //-  Return chemistry-vibrational coupling source term
        //  in the vibrational equation for specie i [kg/m/s3]
        virtual tmp<volScalarField> Scv(const label i) const;

        //-  Return energy removed from free-electrons during impact
        //  ionisation reactions for specie i [kg/m/s3]
        virtual tmp<volScalarField> Siir() const;
        
        //-  Return energy removed from free-electrons during impact
        //  ionisation reactions [kg/m/s3]
        virtual tmp<volScalarField> Siir(const label i) const;


    // I-O

        //- Update properties from given dictionary
        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "hTC2ModelI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
