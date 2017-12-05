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
    selectIdPairsI.H
    selectIdPairs.C
    selectIdPairsIO.C

\*---------------------------------------------------------------------------*/

#ifndef selectIdPairs_H
#define selectIdPairs_H

#include "polyMesh.H"
#include "potential.H"
#include "dictionary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class selectIdPairs Declaration
\*---------------------------------------------------------------------------*/

class selectIdPairs
{
    // Private data
        
        
        // lennard jones
        List<List<scalar> > epsilon_;
        List<List<scalar> > sigma_;

        // buckingham potential
        List<List<scalar> > alpha_;
        List<List<scalar> > rho_;
        List<List<scalar> > C_;
        
        
        bool lennardJones_;
        bool buckinghamPotential_;        
        
        
        label nIds_;




    // Private Member Functions

public:

    // Static data members

    // Constructors

        //- Construct null
        selectIdPairs();

        //- Construct from xs and ys
        selectIdPairs
        (
            const polyMesh& mesh,
            const potential& pot
        );


    // Selectors

    // Destructor

        ~selectIdPairs();

    // Member Functions

            void setConfiguration
            (
                const polyMesh& mesh,
                const potential& pot
            );

        // Access

            inline const List<List<scalar> >& epsilon() const;
            inline const List<List<scalar> >& sigma() const;

            inline const List<List<scalar> >& alpha() const;
            inline const List<List<scalar> >& rho() const;
            inline const List<List<scalar> >& C() const;            
            
            inline const label& nIds() const;

            inline const bool& lennardJones() const;
            inline const bool& buckingham() const;


        // Check

        // Edit

        // Write


    // Member Operators

//         void operator=(const selectIdPairs&);

    // Friend Functions

    // Friend Operators

    // IOstream Operators

//         friend Ostream& operator<<(Ostream&, const selectIdPairs&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "selectIdPairsI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
