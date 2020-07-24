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
    Foam::diffusionModel

Description
    Diffusion model for binary diffusion coefficients

SourceFiles
    diffusionModel.C

\*---------------------------------------------------------------------------*/

#ifndef diffusionModel_H
#define diffusionModel_H

#include "binaryDiffusionModel.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class diffusionModel Declaration
\*---------------------------------------------------------------------------*/

class diffusionModel
{

protected:

    // Protected data

        //- Transport properties dictionary
        const IOdictionary dictTransport_;

        //- Thermo properties dictionary
        const IOdictionary dictThermo_;

        //- Pressure [Pa]
        const volScalarField& p_;

        //- Electron pressure [Pa]
        const volScalarField& pe_;

        //- Temperature [K]
        const volScalarField& T_;

        //- List of species
        wordList species_;

        //- Binary diffusion models
        PtrList<binaryDiffusionModel> DijModels_;

        //- Binary diffusion coefficients
        PtrList<volScalarField> Dij_;


    // Private Member Functions

        //- Disallow copy construct
        diffusionModel(const diffusionModel&);

        //- Disallow default bitwise assignment
        void operator=(const diffusionModel&);


public:

    //- Runtime type information
    TypeName("diffusionModel");

    // Constructors

        //- Construct from components
        diffusionModel
        (
            const word transportPropertiesDictName,
            const word thermoSpeciesDictName,
            const volScalarField& p,
            const volScalarField& pe,
            const volScalarField& T,
            const wordList& species
        );


    // Destructor

        virtual ~diffusionModel()
        {}


    // Member Functions

        inline const volScalarField& D(const label i, const label j) const
        {
            const label iStar = min(i,j);
            const label jStar = max(i,j);
            const label k = species_.size()*iStar+jStar-0.5*iStar*(iStar+1);
            return Dij_[k];
        }

        //- Update the diffusion coefficients
        virtual void update();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
