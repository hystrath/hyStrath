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

#ifndef faceInfo_H
#define faceInfo_H


#include "label.H"
#include "labelField.H"
#include "pTraits.H"

#include "pointInfo.H"
#include "pointInfoList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


class polyMesh;

/*---------------------------------------------------------------------------*\
                           Class faceInfo Declaration
\*---------------------------------------------------------------------------*/

class faceInfo
// :
//         public List< label >
{
    // Private data

        label face_;

//         labelList connectedPoints_;

//         labelList connectedCells_;

        bool interfaceFace_;

        bool acceptedFace_;

        bool zoneFace_;

        bool boundaryFace_;

        bool cyclicFace_;

        bool procFace_;

    // Private Member Functions


        inline void setFace(const label& face);

//         inline void setConnectedPoints(const labelList& points);

//         inline void setConnectedCells
//         (
//             const polyMesh& mesh
//         );


public:

    // Constructors

        //- Construct null
        inline faceInfo();


        //- Construct from mesh and face
        inline faceInfo
        (
//             const polyMesh& mesh,
            const label& face
        );


        //- Construct from components
//         inline faceInfo
//         (
//             const label& face,
//             const labelList& points
//         );


    // Destructor

        ~faceInfo();


    // Member Functions

        inline void setFaceInfo
        (
            const label& face
//             const labelList& points
        );

//         inline void setFaceInfo
//         (
//             const polyMesh& mesh,
//             const label& face
//         );


//         inline void updateFaceInfo
//         (
//             List<pointInfo>& pointData
//         );


//         inline void checkAcceptance
//         (
//             List<pointInfo>& pointData
//         );


        inline void updateFaceBoundaryInfo
        (
            const pointInfoList& pointData,
            const labelList& pointList
        );


        inline void faceOnInterface();

        inline void faceInZone();

        inline void faceOutOfZone();

        inline void acceptFace();

        inline void resetAcceptedFace();

        inline void faceOnBoundary();

        inline void faceOnProcBoundary();

        inline void faceOnCyclicBoundary();



        // Access

            inline const label& faceLabel() const;

            inline const labelList& connectedPoints() const;

            inline const labelList& connectedCells(const polyMesh& mesh);

            inline const bool& isInterfaceFace() const;

            inline const bool& isZoneFace() const;

            inline const bool& isAcceptedFace() const;

            inline const bool& isBoundaryFace() const;

            inline const bool& isProcBoundaryFace() const;

            inline const bool& isCyclicBoundaryFace() const;

        // Check


        // Edit


        // Write

    // Member Operators

    // IOstream Operators

        friend Ostream& operator<<(Ostream&, const faceInfo&);

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "faceInfoI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
