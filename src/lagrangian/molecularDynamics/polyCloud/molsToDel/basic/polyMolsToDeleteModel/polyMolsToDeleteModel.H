/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2007 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Class
    polyMolsToDeleteModel

Description


SourceFiles
    polyMolsToDeleteModelI.H
    polyMolsToDeleteModel.C
    polyMolsToDeleteModelIO.C

\*---------------------------------------------------------------------------*/

#ifndef polyMolsToDeleteModel_H
#define polyMolsToDeleteModel_H

#include "IOdictionary.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "vector.H"
#include "volFields.H"
#include "polyMolecule.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class polyMolsToDeleteModel Declaration
\*---------------------------------------------------------------------------*/

class polyMolsToDeleteModel
{

protected:

    // Protected data

        const fvMesh& mesh_;

        polyMoleculeCloud& molCloud_;


    // Private Member Functions

        void deleteMolFromMoleculeCloud
        (
            polyMolecule& mol
        );

public:

    //- Runtime type information
    TypeName("polyMolsToDeleteModel");

    // Declare runtime constructor selection table
        declareRunTimeSelectionTable
        (
            autoPtr,
            polyMolsToDeleteModel,
            dictionary,
            (
                polyMoleculeCloud& molCloud,
                const dictionary& dict
            ),
            (molCloud, dict)
        );
   
    // Constructors

        //- Construct from components
        polyMolsToDeleteModel
        (
            polyMoleculeCloud& molCloud,
            const dictionary& dict
        );

    // Selectors

        static autoPtr<polyMolsToDeleteModel> New
        (
            polyMoleculeCloud& molCloud,
            const dictionary& dict
        );

    // Destructor

        virtual ~polyMolsToDeleteModel();

    // Member Functions

        // Access
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
