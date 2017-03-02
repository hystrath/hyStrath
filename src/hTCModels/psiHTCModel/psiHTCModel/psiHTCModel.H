/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
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
    Foam::psiHTCModel

Description
    HTC models for compressibility-based thermodynamics

SourceFiles
    psiHTCModelI.H
    psiHTCModel.C
    psiHTCModelNew.C

\*---------------------------------------------------------------------------*/

#ifndef psiHTCModel_H
#define psiHTCModel_H

#include "hTCModel.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "psiReactionThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace hTCModels
{

/*---------------------------------------------------------------------------*\
                     class psiHTCModel Declaration
\*---------------------------------------------------------------------------*/

class psiHTCModel
:
    public hTCModel
{

    // Private Member Functions

        //- Construct as copy (not implemented)
        psiHTCModel(const psiHTCModel&);

        //- Disallow default bitwise assignment
        void operator=(const psiHTCModel&);


public:

    //- Runtime type information
    TypeName("psiHTCModel");


    //- Declare run-time constructor selection tables
    declareRunTimeSelectionTable
    (
        autoPtr,
        psiHTCModel,
        dictionary,
        (
            const word& modelType,
            const fvMesh& mesh
        ),
        (modelType, mesh)
    );


    // Constructors

        //- Construct from components and thermo
        psiHTCModel(const word& modelType, const fvMesh& mesh);



    //- Selector
    static autoPtr<psiHTCModel> New(const fvMesh& mesh);


    //- Destructor
    virtual ~psiHTCModel();


    // Member Functions

        //- Return access to the thermo package
        virtual psiReactionThermo& thermo() = 0;

        //- Return const access to the thermo package
        virtual const psiReactionThermo& thermo() const = 0;

        //- Return tmp of rho
        virtual tmp<volScalarField> rho() const = 0;


    // I-O

        //- Update properties from given dictionary
        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace hTCModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
