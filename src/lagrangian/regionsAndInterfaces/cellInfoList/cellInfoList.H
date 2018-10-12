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


\*---------------------------------------------------------------------------*/

#ifndef cellInfoList_H
#define cellInfoList_H

#include "scalar.H"
#include "vector.H"
#include "point.H"
#include "label.H"
#include "labelField.H"

#include "cellInform.H"
#include "faceInfoList.H"
#include "processorPolyPatch.H"
#include "cyclicPolyPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


class polyMesh;

/*---------------------------------------------------------------------------*\
                           Class cellInfoList Declaration
\*---------------------------------------------------------------------------*/

class cellInfoList
:
        public List< cellInform >
{
    // Private data

    // Private Member Functions

        //- Disallow default bitwise assignment
//         void operator=(const cellInfoList&);

        void setCells();

        void setBoundaryCells( const polyMesh&, const faceInfoList& );

public:

    // Constructors

        //- Construct null
        cellInfoList();

        //- Constructor from list
        cellInfoList( const List<cellInform>& );

        //- Constructor from components
        cellInfoList( const polyMesh&, const faceInfoList&);

        //- Construct as copy
//         inline cellInfoList(const cellInfoList&);


    // Destructor

        ~cellInfoList();

    // Member Functions


        // Access

           void cellsInZone(const labelList&);

//            void cellsInPrefZone(const labelList&);

           void cellsOutOfZone(const labelList&);

//            void cellsOutOfPrefZone(const labelList&);

           void resetAcceptedCells(const labelList&);

        // Check



        // Edit


        // Write

    // Member Operators

        //Note: Used to determine whether to call update.
//         inline bool operator==(const cellInfoList&) const;
//
//         inline bool operator!=(const cellInfoList&) const;

    // IOstream Operators

//         friend Ostream& operator<<(Ostream&, const cellInfoList&);
//         friend Istream& operator>>(Istream&, cellInfoList&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "cellInfoListI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
