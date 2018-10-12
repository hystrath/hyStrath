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
    region

Description
    


SourceFiles


\*---------------------------------------------------------------------------*/

#ifndef region_H
#define region_H

#include "List.H"
#include "scalar.H"
#include "vector.H"
#include "labelField.H"
#include "point.H"

#include "pointInfo.H"
#include "faceInfo.H"
#include "cellInform.H"

#include "pointInfoList.H"
#include "faceInfoList.H"
#include "cellInfoList.H"

#include "interface.H"
#include "processorPolyPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class polyMesh;

/*---------------------------------------------------------------------------*\
                           Class region Declaration
\*---------------------------------------------------------------------------*/

class region
{

    // Private data

        scalar Rsqr_;

        labelList interfacePoints_; //old (points)

        DynamicList<label> pointsToTraverse_;

        DynamicList<label> traversePoints_;

        labelList regionCells_;

        labelList interfaceFaces_; //new (faces)

        word interfaceName_;

        word regionName_;

        List< List<bool> > processorTraverse_;

        List< List< DynamicList<point> > > originPoints_;

        List< List< DynamicList<label> > > procNos_;

        List< DynamicList<label> > originPointpatches_;

        List< List<label> > correspPatchPoints_;

//         DynamicList< DynamicList <label> > pointsToVisualise_;

        labelList acceptedPoints_;

        bool contProcTraverse_;

        List<bool> procTraverse_;


    // Private Member Functions

        inline void traverseZone
        (
            const polyMesh& mesh,
            pointInfoList& pointData,
            const label& originPoint,
            const label& originPointIndex
        );

        inline void traverse
        (
            const polyMesh& mesh,
            pointInfoList& pointData,
            const label& originPointLabel,
            const label& traverseNode,
            const label& originPointIndex
        );

        inline void acceptPointForTraversing
        (
            const label&
        );

        inline bool continueTraverse();

        inline bool contProcTraverse();

        inline void handleBoundaryPoint
        (
            const polyMesh& mesh,
            const pointInfoList& pointData,
            const label& bPoint,
            const point& originPoint,
            const label& originPointLabel,
            const label& originPointIndex,
            const label& procNo
        );

        inline void handleProcBoundaries
        (
            const polyMesh& mesh,
            pointInfoList& pointData
        );

        inline void checkAllProcs(const polyMesh&);

        inline void correspPatchPoints
        (
            const polyMesh& mesh,
            const polyPatch& patch,
            List<label>& patchPoints,
            const List<label>& patchPointFaces,
            const List<label>& facePointIndex
        );

        inline const label& correspondingPatchPoint 
        (
            const polyMesh& mesh,
//             const label& oldPatchPoint,
            const label& faceIndex,
            const label& pointFaceIndex,
            const polyPatch& patch
        ) const;

        void setRegionName(const word&);

        void setRsqr(const scalar&);

        void setZonePoints
        (
            const polyMesh&,
            const labelList&,
            pointInfoList&
        );

        void synchIntProcPoints
        (
            const polyMesh& ,
            pointInfoList& 
        );

        void setInterface
        (
            const polyMesh&,
            const labelList&,
            pointInfoList&
        );

        void setNewInterface
        (
            const polyMesh&,
            faceInfoList&,
            cellInfoList&
        );

        void setCorrPatchPoints
        (
            const polyMesh&,
            const pointInfoList&
        );

        void setRegionPoints
        (
            const polyMesh&,
            pointInfoList&
        );

        void handleBoundaries
        (
            const polyMesh& mesh,
            pointInfoList& pointData
        );

        void setAcceptedPoints(pointInfoList&);

        void setRegionCells
        (
            const polyMesh&,
            pointInfoList&,
            cellInfoList&
        );




public:


    // Constructors

        //- construct null
        region();

        //- construct from layer thickness
        region
        (
            const scalar& R
        );

        //- construct from layer thickness squared and interface points
        region
        (
            const scalar& Rsqr,
            const labelList& interface
        );

    // Destructor

        ~region();

    // Member Functions

        void setRegion
        (
            const polyMesh&,
            const labelList&,
            const labelList&,
            const scalar&,
            const word&,
            pointInfoList&,
            faceInfoList&,
            cellInfoList&
        );

//         inline void setInterface(const labelList& faces);




        // Access

            inline const labelList& acceptedPoints() const;

            inline const labelList& regionCells() const;

            inline const labelList& interfaceFaces() const;

            inline const labelList& interfacePoints() const;

//             inline const DynamicList < DynamicList <label> >& pointsToVisualise() const;

            inline const word& regionName() const;

            inline const word& interfaceName() const;

        // Check


        // Edit


        // Write
        //- Disallow default bitwise assignment
//         void operator=(const region&);


        friend Ostream& operator<<(Ostream&, const region&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "regionI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
