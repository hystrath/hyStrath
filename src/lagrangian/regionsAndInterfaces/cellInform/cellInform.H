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

#ifndef cellInform_H
#define cellInform_H


#include "label.H"
#include "labelField.H"
#include "faceInfo.H"
#include "faceInfoList.H"
#include "pTraits.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


class polyMesh;

/*---------------------------------------------------------------------------*\
                           Class cellInform Declaration
\*---------------------------------------------------------------------------*/

class cellInform
// :
//         public List< faceInfo >
{

    // Private data

        label cell_;

//         labelList connectedFaces_;

        bool acceptedCell_;

        bool prefZoneCell_;

        bool zoneCell_;

        bool boundaryCell_;


    // Private Member Functions

        inline void setCell(const label& cell);

//         inline void setConnectedFaces(const labelList& faces);


public:

    // Constructors

        //- Construct null
        inline cellInform();


        //- Construct from mesh and cell
//         inline cellInform
//         (
// //             const polyMesh& mesh,
//             const label& cell
//         );


        //- Construct from components
        inline cellInform
        (
            const label& cell
//             const labelList& faces
        );


    // Destructor

        ~cellInform();

    // Member Functions

        inline void setCellInfo
        (
            const label& cell
//             const labelList& faces
        );

//         inline void setCellInfo
//         (
//             const polyMesh& mesh,
//             const label& cell
//         );

//         inline void updateCellInfo();


//         inline void checkAcceptance
//         (
//             List<faceInfo>& faceData
//         );


        inline void checkCellBoundary
        (
            const faceInfoList& faceData,
            const labelList& faces
        );

//         inline void checkBoundary
//         (
//             const faceInfoList& faceData
//         );

        inline void cellInZone();

        inline void cellOutOfZone();

        inline void cellInPrefZone();

        inline void cellOutOfPrefZone();

        inline void acceptCell();

        inline void resetAcceptedCell();

        inline void cellOnBoundary();


        // Access

            inline const label& cellLabel() const;

//             inline const labelList& connectedFaces() const;

//             inline const labelList& faces() const;
// 
//             inline const labelList& points() const;

            inline const bool& isPrefZoneCell() const;

            inline const bool& isZoneCell() const;

            inline const bool& isAcceptedCell() const;

            inline const bool& isBoundaryCell() const;

        // Check


        // Edit


        // Write

    // Member Operators

    // IOstream Operators

        friend Ostream& operator<<(Ostream&, const cellInform&);

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "cellInformI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
