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
    Foam::VVModel

Description
    V-V relaxation model for V-V inter-species time relaxation models

SourceFiles
    VVModel.C

\*---------------------------------------------------------------------------*/

#ifndef VVModel_H
#define VVModel_H

#include "VVRelaxationModel.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class VVModel Declaration
\*---------------------------------------------------------------------------*/

class VVModel
{

protected:

    // Protected data

        //- Two-temperature model dictionary
        const IOdictionary dict1_;

        //- Thermophysical properties dictionary
        const IOdictionary dict2_;

        //- Species
        const wordList& species_;

        //- Pressure [Pa]
        const volScalarField& p_;

        //- Translational temperature [K]
        const volScalarField& T_;

        //- Vibrational temperature [K]
        const PtrList<volScalarField>& Tv_;

        //- Number density [m^-3]
        const PtrList<volScalarField>& nD_;


        //- V-V relaxation models
        PtrList<VVRelaxationModel> tauVVijModels_;

        //- Interspecies V-V relaxation time
        PtrList<volScalarField> tauVVij_;


    // Private Member Functions

        //- Disallow copy construct
        VVModel(const VVModel&);

        //- Disallow default bitwise assignment
        void operator=(const VVModel&);


public:

    //- Runtime type information
    TypeName("VVModel");

    // Constructors

        //- Construct from components
        VVModel
        (
            const word& dict1,
            const word& dict2,
            const wordList& species,
            const volScalarField& p,
            const volScalarField& Tt,
            const PtrList<volScalarField>& Tv,
            const PtrList<volScalarField>& nD
        );


    // Destructor

        virtual ~VVModel()
        {}


    // Member Functions

        inline const volScalarField& tauVV(label i, label j) const
        {
            label k = species_.size()*i+j;
            return tauVVij_[k];
        }

        //- Update the diffusion coefficients
        virtual void update();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
