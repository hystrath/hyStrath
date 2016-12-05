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
    Foam::rhoHTCModel

Description
    HTC models for rho-based thermodynamics

SourceFiles
    rhoHTCModelI.H
    rhoHTCModel.C
    rhoHTCModelNew.C

\*---------------------------------------------------------------------------*/

#ifndef rhoHTCModel_H
#define rhoHTCModel_H

#include "hTCModel.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "rhoReactionThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace hTCModels
{

/*---------------------------------------------------------------------------*\
                     class rhoHTCModel Declaration
\*---------------------------------------------------------------------------*/

class rhoHTCModel
:
    public hTCModel
{

    // Private Member Functions

        //- Construct as copy (not implemented)
        rhoHTCModel(const rhoHTCModel&);

        //- Disallow default bitwise assignment
        void operator=(const rhoHTCModel&);


public:

    //- Runtime type information
    TypeName("rhoHTCModel");


    //- Declare run-time constructor selection tables
    declareRunTimeSelectionTable
    (
        autoPtr,
        rhoHTCModel,
        dictionary,
        (
            const word& modelType,
            const fvMesh& mesh
        ),
        (modelType, mesh)
    );


    // Constructors


        //- Construct from components and thermo
        rhoHTCModel(const word& modelType, const fvMesh& mesh);



    //- Selector
    static autoPtr<rhoHTCModel> New(const fvMesh& mesh);


    //- Destructor
    virtual ~rhoHTCModel();


    // Member Functions

        // Access functions

            //- Access hTC dict
             inline const dictionary& coeff() const;


        //- Return access to the thermo package
        virtual rhoReactionThermo& thermo() = 0;

        //- Return const access to the thermo package
        virtual const rhoReactionThermo& thermo() const = 0;

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
