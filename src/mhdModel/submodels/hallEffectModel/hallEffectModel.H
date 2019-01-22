/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2014 OpenFOAM Foundation
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
    Foam::mhd:hallEffectModel

Description
    Model to supply absorption and emission coefficients for radiation
    modelling

\*---------------------------------------------------------------------------*/

#ifndef mhdHallEffectModel_H
#define mhdHallEffectModel_H

#include "IOdictionary.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "volFields.H"
#include "Vector2D.H"
#include "mhdModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace mhd
{

/*---------------------------------------------------------------------------*\
                   Class conductivityModel Declaration
\*---------------------------------------------------------------------------*/

class hallEffect
{

protected:

    // Protected data

        //- mhd model dictionary
        const volVectorField& B_;

        //- Reference to the fvMesh
        const fvMesh& mesh_;


public:

    //- Runtime type information
    TypeName("hallEffectModel");

    //- Declare runtime constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            hallEffectModel,
            mhdModel,
            (
                const volVectorField& B,
                const fvMesh& mesh
            ),
            (dict, mesh)
        );


    // Constructors

        //- Construct from components
        hallEffectModel
        (
            const volVectorField& B,
            const fvMesh& mesh
        );


    //- Selector
    static autoPtr<hallEffectModel> New
    (
        const volVectorField& B,
        const fvMesh& mesh
    );


    //- Destructor
    virtual ~hallEffectModel();


    // Member Functions

        // Access

            //- Reference to the mesh
            inline const fvMesh& mesh() const
            {
                return mesh_;
            }

            //- Reference to the dictionary
            inline const bolVector& V() const
            {
                return B_;
            }



                //- hallParameter
                virtual volTensorField () const = 0;


            //- Correct
            virtual void correct();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace radiation
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

