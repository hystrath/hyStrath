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
    Foam::eVModel

Description
    e-V relaxation model

SourceFiles
    eVModel.C

\*---------------------------------------------------------------------------*/

#ifndef eVModel_H
#define eVModel_H

#include "eVRelaxationModel.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class eVModel Declaration
\*---------------------------------------------------------------------------*/

class eVModel
{

protected:

    // Protected data

        //- Two-temperature model dictionary
        const IOdictionary dict2T_;

        //- Thermophysical properties dictionary
        const IOdictionary dictThermoPhy_;

        //- Species
        const wordList& species_;

        //- Pressure [Pa]
        const volScalarField& p_;

        //- Vibrational temperature [K]
        const volScalarField& Tv_;

        //- e-V relaxation models
        PtrList<eVRelaxationModel> taueViModels_;

        //- e-V relaxation time
        PtrList<volScalarField> taueVi_;


    // Private Member Functions

        //- Disallow copy construct
        eVModel(const eVModel&);

        //- Disallow default bitwise assignment
        void operator=(const eVModel&);


public:

    //- Runtime type information
    TypeName("eVModel");

    // Constructors

        //- Construct from components
        eVModel
        (
            const word& dict2T,
            const word& dictThermoPhy,
            const wordList& species,
            const volScalarField& p,
            const volScalarField& Tv
        );


    // Destructor

        virtual ~eVModel()
        {}


    // Member Functions

        inline const volScalarField& taueV(label i) const
        {
            return taueVi_[i];
        }

        //- Update the e-V relaxation times
        virtual void update();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
