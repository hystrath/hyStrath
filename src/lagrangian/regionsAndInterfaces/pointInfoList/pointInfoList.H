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

#ifndef pointInfoList_H
#define pointInfoList_H

#include "scalar.H"
#include "vector.H"
#include "point.H"
#include "label.H"
#include "labelField.H"

#include "pointInfo.H"
#include "processorPolyPatch.H"
#include "cyclicPolyPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


class polyMesh;

/*---------------------------------------------------------------------------*\
                           Class pointInfoList Declaration
\*---------------------------------------------------------------------------*/

class pointInfoList
:
        public List< pointInfo >
{
    // Private data

//         List< List <label> > patchPoints_;

    // Private Member Functions

        //- Disallow default bitwise assignment
//         void operator=(const pointInfoList&);

        void setPoints( const polyMesh& );

        void setBoundaryPoints( const polyMesh& );

//         inline void findGlobalPoint
//         (
//             const polyMesh& mesh,
//             const polyPatch& patch,
//             const label& localPoint
// //             const label& localPointIndex
//         );

        //- Disallow default bitwise assignment
        void operator=(const pointInfoList&);

        //- Disallow default bitwise copy construct
        pointInfoList(const pointInfoList&);

public:

    // Constructors

        //- Construct null
        pointInfoList();

        //- Constructor from list
        pointInfoList( const List<pointInfo>& );

        //- Constructor from components
        pointInfoList( const polyMesh& );

        //- Construct as copy
//         inline pointInfoList(const pointInfoList&);


    // Destructor

        ~pointInfoList();

    // Member Functions


        // Access

            void clearInterfacePoints();

            void unSwitchProcPoints();

            void unVisitPoints();

            void unAcceptPoints();

            void clearZonePoints();

        // Check

            void switchProcPoints(const labelList&);


        // Edit


        // Write

    // Member Operators

        //Note: Used to determine whether to call update.
//         inline bool operator==(const pointInfoList&) const;
//
//         inline bool operator!=(const pointInfoList&) const;

    // IOstream Operators
//         friend Ostream& operator<< (Ostream&, const pointInfoList&);
// 
//         friend Istream& operator>> (Istream&, pointInfoList&);




//         friend Ostream& operator<<(Ostream&, const pointInfoList&);
//         friend Istream& operator>>(Istream&, pointInfoList&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "pointInfoListI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
