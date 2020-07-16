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
    noExclusionsI.H
    noExclusions.C
    noExclusionsIO.C

\*---------------------------------------------------------------------------*/

#ifndef noExclusions_H
#define noExclusions_H

#include "polyMoleculeCloud.H"
#include "exclusionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class noExclusions Declaration
\*---------------------------------------------------------------------------*/

class noExclusions
:
    public exclusionModel
{

private:

    // Private data

        dictionary propsDict_;




public:

    //- Runtime type information
    TypeName("noExclusions");

    // Static data members

    // Constructors

        //- construct from components
        noExclusions
        (
            const polyMesh& mesh,
            polyMoleculeCloud& molCloud,
            const dictionary& dict
        );

    // Selectors


    // Destructor

        ~noExclusions();

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

//         void operator=(const noExclusions&);

    // Friend Functions

    // Friend Operators

    // IOstream Operators

//         friend Ostream& operator<<(Ostream&, const noExclusions&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
