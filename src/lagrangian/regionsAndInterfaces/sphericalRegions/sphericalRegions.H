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

#ifndef sphericalRegions_H
#define sphericalRegions_H

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


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

 class polyMesh;

/*---------------------------------------------------------------------------*\
                           Class sphericalRegions Declaration
\*---------------------------------------------------------------------------*/

class sphericalRegions
:
        public List< region >
{

    // Private data

        labelList remainingCells_;

        labelList sphericalRegionsCells_;


    // Private Member Functions

        void setRegionsCells
        (
            const labelList&
        );

        void setRemainingCells
        (
            const labelList&
        );


public:

    // Constructors

        //- Construct from components: polyMesh, list of pointInfo,
        //   list of faces (interface), list of cells (zone) and
        //   a list of sub layer thicknesses.


        sphericalRegions
        (
            const polyMesh&,
            pointInfoList&,
            faceInfoList&,
            cellInfoList&,
            const labelList&,
            const labelList&,
            const List<word>&,
            const List<scalar>&
        );


    // Destructor

        ~sphericalRegions();

    // Member Functions


        // Access

            inline const labelList& remainingCells() const;

            inline const labelList& regionsCells() const;

        // Check

        // Edit

        // Write

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "sphericalRegionsI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
