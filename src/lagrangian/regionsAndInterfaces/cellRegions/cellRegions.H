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

#ifndef cellRegions_H
#define cellRegions_H

#include "List.H"
#include "scalar.H"
#include "vector.H"
#include "labelField.H"

#include "pointInfo.H"
#include "faceInfo.H"
#include "cellInform.H"

#include "pointInfoList.H"
#include "faceInfoList.H"
#include "cellInfoList.H"

#include "region.H"
#include "labelIOList.H"
#include "IOobjectList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

 class polyMesh;

/*---------------------------------------------------------------------------*\
                           Class cellRegions Declaration
\*---------------------------------------------------------------------------*/

class cellRegions
:
        public List< region >
{

    // Private data

        scalar Rsqr_;

        //- List of processor cell addressing lists
        labelList cellProcAddressing_;



//         labelList remainingCells_;
// 
//         const labelList& zoneCells_;


    // Private Member Functions

//         void setRegionsCells
//         (
//             const labelList&
//         );
// 
//         void setRemainingCells
//         (
//             const labelList&
//         );


public:

    // Constructors

        //- Construct from components: polyMesh, list of pointInfo,
        //   list of faces (interface), list of cells (zone) and
        //   a list of sub layer thicknesses.


        cellRegions
        (
            const polyMesh&,
            pointInfoList&,
            faceInfoList&,
            cellInfoList&,
            const labelList&,
            const scalar&
        );


    // Destructor

        ~cellRegions();

    // Member Functions


        // Access
/*
            inline const labelList& remainingCells() const;

            inline const labelList& regionsCells() const;*/

        // Check

        // Edit

        // Write

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "cellRegionsI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
