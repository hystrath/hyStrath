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


Description

SourceFiles
    selectIdsI.H
    selectIds.C
    selectIdsIO.C

\*---------------------------------------------------------------------------*/

#ifndef selectIds_H
#define selectIds_H

#include "polyMesh.H"
#include "constantMoleculeProperties.H"
#include "dictionary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class selectIds Declaration
\*---------------------------------------------------------------------------*/

class selectIds
{
    // Private data

        List<label> molIds_;

        List<word> molIdNames_;

    // Private Member Functions

public:

    // Static data members

    // Constructors

        //- Construct null
        selectIds();

        selectIds
        (
            const constantMoleculeProperties& cP,
            const dictionary& dict
        );

        selectIds
        (
            const constantMoleculeProperties& cP,
            const dictionary& dict,
            const word& molIdsHeader
        );

    // Selectors

    // Destructor

        ~selectIds();

    // Member Functions

        // Access

            inline List<label> molIds();
            inline List<word> molIdNames();
//             inline labelList& molIds();


        // Check

        // Edit

        // Write


    // Member Operators

//         void operator=(const selectIds&);

    // Friend Functions

    // Friend Operators

    // IOstream Operators

//         friend Ostream& operator<<(Ostream&, const selectIds&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "selectIdsI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
