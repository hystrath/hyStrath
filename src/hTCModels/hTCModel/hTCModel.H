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
    Foam::hTCModel

Description
    Base class for hTC models

SourceFiles
    hTCModel.C

\*---------------------------------------------------------------------------*/

#ifndef hTCModel_H
#define hTCModel_H

#include "IOdictionary.H"
#include "turbulenceModel2.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class hTCModel Declaration
\*---------------------------------------------------------------------------*/

class hTCModel
:
    public IOdictionary
{

protected:

    // Protected data

        //- Reference to the turbulence model
        compressible::turbulenceModel2* turbulencePtr_;

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
        hTCModel(const hTCModel&);

        //- Disallow default bitwise assignment
        void operator=(const hTCModel&);


public:

    //- Runtime type information
    TypeName("hTCModel");


    // Constructors


        //- Construct from components
        hTCModel(const word& modelType, const fvMesh& mesh);


    //- Destructor
    virtual ~hTCModel();


    // Member Functions

        // Access


            //- Return const access to the mesh database
            inline const fvMesh& mesh() const;

            //- Return const access to phi
            inline const surfaceScalarField& phi() const;

            //- Return const access to rho
            virtual tmp<volScalarField> rho() const = 0;

            //- Return access to turbulence
            inline const compressible::turbulenceModel2& turbulence() const;

            //- Set turbulence
            inline void setTurbulence
            (
                compressible::turbulenceModel2& turbModel
            );

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


    // I-O

        //- Update properties from given dictionary
        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "hTCModelI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
