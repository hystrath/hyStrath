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
    cyclicBoundary

Description
    This is the definition of a cyclic boundary - data is read in from the 
    cyclicBoundaries dictionary (see constructCyclicBoundaries class and 
    createCyclicBoundaries utility for more info). 
    
    Only linear translations has been implemented between coupled cyclis,
    but it is possible to introduce rotational translations as well in the future. 
    This was present in 1.7.x but it was discontinued until it is required.


SourceFiles
    cyclicBoundaryI.H
    cyclicBoundary.C
    cyclicBoundaryIO.C

\*---------------------------------------------------------------------------*/

#ifndef cyclicBoundary_H
#define cyclicBoundary_H

#include "IOdictionary.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "vector.H"
#include "volFields.H"
#include "polyMesh.H"
#include "Time.H"

#include "writeTimeData.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class cyclicBoundary Declaration
\*---------------------------------------------------------------------------*/

class cyclicBoundary
{

protected:

    // Protected data

        const fvMesh& mesh_;

        //- subDictionary containing the properties
        dictionary boundaryDict_;

        //- name of this patch
        word patchName_;
        
        //- id on boundaryMesh
        label patchId_;

        // normal vector (pointing outwards)
        vector nF_;
        
        // centroid of patch 
        vector centroid_;
        
        // 4 points (assuming a rectangular patch)
        List<vector> boundaryPoints_;

        // total number of faces on patch
        label nFaces_;
        
        // separationVector direction from this patch to neighbour
        vector separationVector_;        

        //- name of neighbour patch
        word patchNameN_;        
        label patchNId_;

        // angle between coupled patches
        scalar theta_;
        bool rotate_;
// 
//         vector rotationStartPoint_;
//         vector rotationEndPoint_;
//         vector rotationAxis_;
//         vector rotationPt_;
// 
//         // rotational tensors
//         tensor RAB_;
//         tensor RBA_;

        // boundary normals
//         vector nA_;
//         vector nB_;

        // boundary mid-points
//         vector pA_;
//         vector pB_;

//         label nFaces_;

        //- faces on the cyclicBoundary
//         labelList faces_;

//         labelList coupledFacesA_;
//         labelList coupledFacesB_;

        //- cells corresponding to the cyclicBoundary faces
//         labelList cells_;

//         labelList cellsA_;
//         labelList cellsB_;

    // Private Member Functions

//         void getCoupledFaces(const polyMesh& mesh);

//         pointField points
//         (
//             const List<vector>& points,
//             const vector& normal,
//             const scalar& offset
//         );
//         void readInfo();

public:

    //- Runtime type information
    TypeName("cyclicBoundary");

    // Declare runtime constructor selection table
        declareRunTimeSelectionTable
        (
            autoPtr,
            cyclicBoundary,
            dictionary,
            (
                Time& t,
                const polyMesh& mesh,
                const dictionary& dict
            ),
            (t, mesh, dict)
        );
   
    // Constructors

        //- Construct from components
        cyclicBoundary
        (
            Time& t,
            const polyMesh& mesh,
            const dictionary& dict
        );


    // Selectors

        static autoPtr<cyclicBoundary> New
        (
            Time& t,
            const polyMesh& mesh,
            const dictionary& dict
        );


    // Destructor

        virtual ~cyclicBoundary();


    // Member Functions
    
            
        bool isPatchNeighbour(const word& patchNameNeighbour);
            
        void setCoupledPatchInfo
        (
            const vector& centroid,
            const vector& nF
        );

        // Access

            //- return the name of the control patch
            const word& patchName() const;
            const label& patchId() const;
            const vector& separationVector() const;            
            const List<vector>& boundaryPoints() const;
            const vector& normal() const;             
            const vector& centroid() const;               

            // neighbour
            const word& patchNameN() const;
            const label& patchIdN() const;           
            
            // parallel reconstructed number of faces
//             const label& nFaces() const;
// 
//             //- retrun faces on patch
//             const labelList& controlPatch() const;
//             const labelList& controlPatchA() const;
//             const labelList& controlPatchB() const;

//             const vector& normalA() const;
//             const vector& normalB() const;

//             const vector& translation() const;

            //- return the adjacent cells
//             const labelList& controlZone() const;
//             const labelList& controlZoneA() const;
//             const labelList& controlZoneB() const;





            const scalar& theta() const;
            const bool& rotate() const;
//             const vector& rotationAxis() const;
//             const vector& rotationAxisStartPt() const;
// 
//             const vector& rotationPt() const;
//             const tensor& RAB() const;
//             const tensor& RBA() const;

    
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
