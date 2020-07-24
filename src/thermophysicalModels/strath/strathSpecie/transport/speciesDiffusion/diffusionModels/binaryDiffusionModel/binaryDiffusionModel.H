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
    Foam::binaryDiffusionModel

Description
    An abstract base class for binary diffusion models.

SourceFiles
    binaryDiffusionModel.C

\*---------------------------------------------------------------------------*/

#ifndef binaryDiffusionModel_H
#define binaryDiffusionModel_H

#include "specie.H"
#include "IOdictionary.H"
#include "typeInfo.H"
#include "runTimeSelectionTables.H"
#include "scalarField.H"
#include "volFieldsFwd.H"
#include "surfaceFieldsFwd.H"
#include "dimensionedScalar.H"
#include "tmp.H"
#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class binaryDiffusionModel Declaration
\*---------------------------------------------------------------------------*/

class binaryDiffusionModel
{

protected:

    // Protected data

        //- Species name
        const word& name1_, name2_;

        //- Thermo properties dictionary
        const dictionary& dictThermo_;

        //- Transport properties dictionary
        const dictionary& dictTransport_;

        //- Pressure [Pa]
        const volScalarField& p_;

        //- Electron pressure [Pa]
        const volScalarField& pe_;

        //- Temperature [K]
        const volScalarField& T_;

        //- Integer to give the type of collision interaction
        //  0: 2 neutrals, 1: 1 neutral-ion/e-, 2: 2 ions, 3: ion-electron
        //  4: 2 electrons
        label collisionType_;


    // Private Member Functions

        //- Disallow copy construct
        binaryDiffusionModel(const binaryDiffusionModel&);

        //- Disallow default bitwise assignment
        void operator=(const binaryDiffusionModel&);


public:

    //- Runtime type information
    TypeName("binaryDiffusionModel");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            binaryDiffusionModel,
            dictionary,
            (
                const word& name1,
                const word& name2,
                const dictionary& dictThermo,
                const dictionary& dictTransport,
                const volScalarField& p,
                const volScalarField& pe,
                const volScalarField& T
            ),
            (name1, name2, dictThermo, dictTransport, p, pe, T)
        );


    // Selectors

        //- Return a reference to the selected diffusion model
        static autoPtr<binaryDiffusionModel> New
        (
            const word& name1,
            const word& name2,
            const dictionary& dictThermo,
            const dictionary& dictTransport,
            const volScalarField& p,
            const volScalarField& pe,
            const volScalarField& T
        );


    // Constructors

        //- Construct from components
        binaryDiffusionModel
        (
            const word& name1,
            const word& name2,
            const dictionary& dictThermo,
            const dictionary& dictTransport,
            const volScalarField& p,
            const volScalarField& pe,
            const volScalarField& T
        );


    // Destructor

        virtual ~binaryDiffusionModel()
        {}


    // Member Functions

        //- Binary diffusion coefficient [m^2/s]
        virtual tmp<volScalarField> D() const;

        //- Binary diffusion coefficient for patch [m^2/s]
        //  with addition of the electron pressure
        virtual tmp<scalarField> D
        (
            const scalarField& p,
            const scalarField& pe,
            const scalarField& T,
            const label patchi
        ) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
