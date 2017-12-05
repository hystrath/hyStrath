/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2005 OpenCFD Ltd.
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


Description

SourceFiles
    excludeAllI.H
    excludeAll.C
    excludeAllIO.C

\*---------------------------------------------------------------------------*/

#ifndef excludeAll_H
#define excludeAll_H

#include "polyMoleculeCloud.H"
#include "exclusionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class excludeAll Declaration
\*---------------------------------------------------------------------------*/

class excludeAll
:
    public exclusionModel
{

private:

    // Private data

        dictionary propsDict_;




public:

    //- Runtime type information
    TypeName("excludeAll");

    // Static data members

    // Constructors

        //- construct from components
        excludeAll
        (
            const polyMesh& mesh,
            polyMoleculeCloud& molCloud,
            const dictionary& dict
        );

    // Selectors


    // Destructor

        ~excludeAll();

    // Member Functions

        // Access

        // Check

        // Edit
        void initialiseExclusions();

        bool excludeMolecules
        (
            polyMolecule* molI,
            polyMolecule* molJ
        );
        
        bool excludeSites
        (
            polyMolecule* molI,
            polyMolecule* molJ,
            const label& siteI,
            const label& siteJ
        );
        

        // Write





    // Member Operators

//         void operator=(const excludeAll&);

    // Friend Functions

    // Friend Operators

    // IOstream Operators

//         friend Ostream& operator<<(Ostream&, const excludeAll&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
