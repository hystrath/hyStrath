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
    molConfig

Description

SourceFiles
    molConfigI.H
    molConfig.C
    molConfigIO.C

\*---------------------------------------------------------------------------*/

#ifndef faceInfoList_H
#define faceInfoList_H

#include "scalar.H"
#include "vector.H"
#include "point.H"
#include "label.H"
#include "labelField.H"

#include "faceInfo.H"
#include "processorPolyPatch.H"
#include "cyclicPolyPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


class polyMesh;

/*---------------------------------------------------------------------------*\
                           Class faceInfoList Declaration
\*---------------------------------------------------------------------------*/

class faceInfoList
:
        public List< faceInfo >
{
    // Private data

    // Private Member Functions

        //- Disallow default bitwise assignment
//         void operator=(const faceInfoList&);

        void setFaces();

        void setBoundaryFaces( const polyMesh& );

public:

    // Constructors

        //- Construct null
        faceInfoList();

        //- Constructor from list
        faceInfoList( const List<faceInfo>& );

        //- Constructor from components
        faceInfoList( const polyMesh& );

        //- Construct as copy
//         inline faceInfoList(const faceInfoList&);


    // Destructor

        ~faceInfoList();

    // Member Functions


        // Access

            void facesOnInterface( const labelList& );

            void removeAcceptedFaces();
        // Check



        // Edit


        // Write

    // Member Operators

        //Note: Used to determine whether to call update.
//         inline bool operator==(const faceInfoList&) const;
//
//         inline bool operator!=(const faceInfoList&) const;

    // IOstream Operators

//         friend Ostream& operator<<(Ostream&, const faceInfoList&);
//         friend Istream& operator>>(Istream&, faceInfoList&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "faceInfoListI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
