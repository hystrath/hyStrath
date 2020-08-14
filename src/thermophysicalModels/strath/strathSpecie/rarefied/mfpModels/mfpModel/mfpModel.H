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
    Foam::mfpModel

Description
    Mean free path model for mfp abstract class

SourceFiles
    mfpModel.C

\*---------------------------------------------------------------------------*/

#ifndef mfpModel_H
#define mfpModel_H

#include "basicMfpModel.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class mfpModel Declaration
\*---------------------------------------------------------------------------*/

class mfpModel
{

protected:

    // Protected data

        //- Transport properties dictionary
        const IOdictionary dict_;

        //- Thermophysical properties dictionary
        const IOdictionary dictThermoPhy_;

        //- Species
        const wordList& species_;

        //- Pressure [Pa]
        const volScalarField& p_;

        //- Heavy particle temperature [K]
        const volScalarField& T_;


        //- Mean free path models
        PtrList<basicMfpModel> mfpModels_;

        //- Mean free path
        PtrList<volScalarField> mfp_;


    // Private Member Functions

        //- Disallow copy construct
        mfpModel(const mfpModel&);

        //- Disallow default bitwise assignment
        void operator=(const mfpModel&);


public:

    //- Runtime type information
    TypeName("mfpModel");

    // Constructors

        //- Construct from components
        mfpModel
        (
            const word transportPropertiesDictName,
            const word thermoSpeciesDictName,
            const wordList& species,
            const volScalarField& p,
            const volScalarField& Tt
        );


    // Destructor

        virtual ~mfpModel()
        {}


    // Member Functions

        //- Returns the constant mfp for species i
        inline const volScalarField& mfp(label i) const
        {
            return mfp_[i];
        }

        //- Returns the mfp for species i
        inline volScalarField& mfp(label i)
        {
            return mfp_[i];
        }

        //- Update the mean free path of each species
        virtual void update();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
