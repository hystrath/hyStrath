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
    cellInteractions

Description



SourceFiles
    cellInteractionsI.H
    cellInteractions.C

\*---------------------------------------------------------------------------*/

#ifndef cellInteractions_H
#define cellInteractions_H

#include "vector.H"
#include "fvMesh.H"
#include "referredCell.H"
#include "patchLayer.H"
#include "cyclicBoundaries.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class cellInteractions Declaration
\*---------------------------------------------------------------------------*/

template<class ParticleType>
class cellInteractions
{
    // Private data

        const polyMesh& mesh_;

        const cyclicBoundaries& cyclics_;

        const reducedUnits& redUnits_;
        
        word fieldName_;

        // dummy cloud
        Cloud<ParticleType> cloud_;

        // the mesh is optimised when the cell is a cube of side-length rCut
        // (or slightly larger than rCut)
        bool optimised_;

        // distance over which interaction lists will be constructed
        scalar rCut_;

        // the optimised direct interaction list (avoids double-count looping)
        // outer list is the size of this mesh
        // inner list are the unique neighbouring cells to interact with
        labelListList dil_;

        // the full interaction list (does not cater of double-count looping like DIL)
        labelListList fil_;


        // outer list is mesh_.nCells()
        // inner list are the ids to those cells in the DIL interaction of cellI (not including cellI)
        // that abide with dil_ structure 
        labelListList inverseDIL_;


        // outer list is processor numbers, 
        // inner list is referred ids for receiving data. 
        // The stored ids correspond to the data member: refCells_
        labelListList recRefIds_;

        // outer list is processor numbers, 
        // inner list is source cell labels for sending data
        // This processor should know which real particles in its source cells 
        // to send to other processors
        labelListList sendSrcCells_;
        
        // list of referred cells 
        List<referredCell> refCells_;

        // outer list is the mesh cells
        // inner list is the ids to the refCells
        // own processor only - i.e. periodic cells
        // if a real-cell is known, can you find which referred cells are
        //  imaged from this source cell
        labelListList sourceCellToRefs_;

        
        // outer list is the mesh cells
        // inner list is the ids to the refCells 
        // full real-to-referred interaction list - from a real cell, obtain the full list
        // of referred cells within interactions (ids only)
        labelListList inverseFRIL_;
        
        // similar to refCells but with molecules per referred cell
        // outer list is referred cells 
        // inner list is referred particles within referred cell
        List<DynamicList<ParticleType*> > refCellsParticles_;

        // each (real) cell on the mesh stores a list of the full referred list of 
        // particles for interaction
        // outer list is mesh_.nCells()
        // inner list is list of referred particles within interaction of real cell, cellI
        List<DynamicList<ParticleType*> > ripl_;
        
        //- container of referred particles in no particular order (local mesh)
        IDLList<ParticleType> referredCloud_;

        label nCoupledPatches_;
        label nCyclicProcessorBoundaries_;
        label nProcessorBoundaries_;
        label nCyclicBoundaries_;

        bool write_;
        bool writeReferredCells_;
        

    // Private Member Functions

            void buildInteractionLists();

            void buildOptimisedDIL();

            void buildComplexDIL();
            
            void buildInverseDirectInteractionList();            
            
            void buildFullInteractionList();
            
            void buildReferredCells();
            

            pointField boundingCellPoints(const label& cellI);

            bool includeCell
            (
                const List<DynamicList<label> >& cells,
                const label& cellI,
                const label& cellN
            );

            bool interactingCells
            (
                const label& cellI,
                const label& cellN,
                const scalar& offset
            );

            boundedBox cellToBoundBox(const label& cellI);

            labelList getProcessorSourceCells
            (
                const polyPatch& patch
            );             
            
            boundedBox faceToBoundBox
            (
                const label& faceI,
                const scalar& smallOffset
            );
           
            labelList returnNeighbouringCells(const label& cellI);
  
            //debug
            void checkForOverlaps();
            
public:

    // Constructors

        //- Construct from mesh 
        cellInteractions
        (
            const polyMesh& mesh,
            const reducedUnits& rU,           
            const cyclicBoundaries& cyclics,
            const scalar& rCut,
            const word& fieldName
        );


        // construct as NULL - dummy inputs
        cellInteractions
        (
            const polyMesh& mesh,
            const reducedUnits& rU,  
            const cyclicBoundaries& cyclics,
            const label& dummy
        );

    // Destructor

        ~cellInteractions();

    // Member Functions
        void checkMesh(const scalar& rCut);

        void addParticle
        (
            ParticleType* particle,
            DynamicList<ParticleType*>& newRefParticles,
            DynamicList<label>& refIds
        );
        
        void deleteParticle
        (
            ParticleType* particle,
            DynamicList<ParticleType*>& newRefParticles,
            DynamicList<label>& refIds
        );
        
        void setReferredParticles
        (
            const List<DynamicList<ParticleType*> >& cellOccupancy
        );              

        void setRIPL();
        
        // Access
        
            inline const labelListList& inverseDIL() const;

            inline const labelListList& dil() const;

            inline const labelListList& fil() const;

            inline const List<referredCell>& refCells() const;
            inline List<referredCell>& refCells();

            inline const List<DynamicList<ParticleType*> >& refCellsParticles() const;
            inline List<DynamicList<ParticleType*> >& refCellsParticles();

            inline const labelListList& sourceCellToRefs() const;

            inline const labelListList& inverseFRIL() const;            
            
            inline const List<DynamicList<ParticleType*> >& ripl() const;

            inline IDLList<ParticleType>& referredCloud();
               
        // Check

        // Edit

        // Write

            bool write();

            void writeReferredCells();

    // Member Operators

    // Friend Functions

    // Friend Operators
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "cellInteractionsI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "cellInteractions.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
