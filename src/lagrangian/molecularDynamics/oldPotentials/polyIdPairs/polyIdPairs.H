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
    polyIdPairsI.H
    polyIdPairs.C
    polyIdPairsIO.C

\*---------------------------------------------------------------------------*/

#ifndef polyIdPairs_H
#define polyIdPairs_H

#include <iostream>
#include <vector>

#include "polyMesh.H"
#include "potential.H"
#include "dictionary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class polyIdPairs Declaration
\*---------------------------------------------------------------------------*/

class polyIdPairs
{
    // Private data
        
        List<List<List<scalar> > > coeffVals_;
        List<word> coeffNames_;
        List<label> coeffNumIds_;
        label nIds_;
        label coeffSize_;
        std::string coeffType_;


    // Private Member Functions

public:

    // Static data members

    // Constructors

        //- Construct null
        polyIdPairs();

        //- Construct from xs and ys
        polyIdPairs
        (
            const polyMesh& mesh,
            const potential& pot
        );


    // Selectors

    // Destructor

        ~polyIdPairs();

    // Member Functions

        // Access

            inline const label& nIds() const;
            inline const label& coeffSize() const;
            inline const List<List<List<scalar> > >& coeffVals() const;
            inline const List<word>& coeffNames() const;
            inline const List<label>& coeffNumIds() const;
            inline const std::string& coeffType() const;


        // Check

        // Edit

        // Write


    // Member Operators

//         void operator=(const polyIdPairs&);

    // Friend Functions

    // Friend Operators

    // IOstream Operators

//         friend Ostream& operator<<(Ostream&, const polyIdPairs&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "polyIdPairsI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
