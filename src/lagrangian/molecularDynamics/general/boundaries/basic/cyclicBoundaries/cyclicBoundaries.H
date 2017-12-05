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
    cyclicBoundaries

Description




\*---------------------------------------------------------------------------*/

#ifndef cyclicBoundaries_H
#define cyclicBoundaries_H

#include "List.H"
#include "cyclicBoundary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class cyclicBoundaries Declaration
\*---------------------------------------------------------------------------*/

class cyclicBoundaries
{
    // Private data

        Time& time_;
        
        //- The entire dictionary (containing multiple subDictionaries)
        IOdictionary boundariesDict_;
        
        //- cyclic boundary models
        PtrList<entry> cyclicBoundaryList_;
        List<word> cyclicBoundaryNames_;
        List<label> cyclicBoundaryIds_;
        List<fileName> cMFixedPathNames_;
        List< autoPtr<cyclicBoundary> > cyclicBoundaryModels_;
        List<label> cyclicBoundaryToModelId_;


public:

    // Constructors
    
        //- Null Constructor 
        cyclicBoundaries
        (
            Time& t,
            const polyMesh& mesh
        );

        //- Constructor 
        cyclicBoundaries
        (
            Time& t,
            const polyMesh& mesh,
            const label& dummy
        );

        ~cyclicBoundaries();


    // Member Functions

//         const label& nCyclicBoundaryModels() const;


        inline List< autoPtr<cyclicBoundary> >& cyclicBoundaryModels();
        inline const List< autoPtr<cyclicBoundary> >& cyclicBoundaryModels() const;
        inline const List<word>& cyclicBoundaryNames() const;
        inline const List<label>& cyclicBoundaryIds() const;
        inline const List<label>& cyclicBoundaryToModelIds() const;



    // Ostream Operator

//         friend Ostream& operator<<(Ostream&, const cyclicBoundaries&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "cyclicBoundariesI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
