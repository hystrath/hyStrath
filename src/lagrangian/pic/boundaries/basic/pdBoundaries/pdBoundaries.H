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
    pdBoundaries

Description

	Stores all the information for the boundary models. Models are defined
    within the boundariesDict. This class selects & builds the models 
    automatically. We assume that boundaries are divided into three:

    1. Patch models -- actions on instantaneous collisions with a patch boundary 

    2. cyclic boundary models -- actions on instantaneous collisions with 
                                a cyclic boundary 

    3. General boundary models -- time-dependent actions on molecules close to
                                  the boundary


\*---------------------------------------------------------------------------*/

#ifndef pdBoundaries_H
#define pdBoundaries_H

#include "List.H"

#include "pdCyclicBoundary.H"
#include "pdPatchBoundary.H"
#include "pdGeneralBoundary.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class pdBoundaries Declaration
\*---------------------------------------------------------------------------*/

class pdBoundaries
{
    // Private data


        Time& time_;

        //- The entire dictionary (containing multiple subDictionaries)
        IOdictionary pdBoundariesDict_;

        //- number of boundary models
        label nPatchBoundaryModels_;
        label nCyclicBoundaryModels_;
        label nGeneralBoundaryModels_;

        //- patch models
        PtrList<entry> patchBoundaryList_;
        List<word> patchBoundaryNames_;
        List<label> patchBoundaryIds_;
        List<fileName> pBFixedPathNames_;
        List< autoPtr<pdPatchBoundary> > patchBoundaryModels_;
        List<label> patchToModelId_;

        //- cyclic boundary models
        PtrList<entry> cyclicBoundaryList_;
        List<word> cyclicBoundaryNames_;
        List<label> cyclicBoundaryIds_;
        List<fileName> cMFixedPathNames_;
        List< autoPtr<pdCyclicBoundary> > cyclicBoundaryModels_;
        List<label> cyclicBoundaryToModelId_;

        //- general boundary models 
        PtrList<entry> generalBoundaryList_;
        List<word> generalBoundaryNames_;
        List<label> generalBoundaryIds_;
        List<fileName> gMFixedPathNames_;
        List< autoPtr<pdGeneralBoundary> > generalBoundaryModels_;





public:

    // Constructors
    
        //- Null Constructor 
        pdBoundaries
        (
            Time& t,
            const polyMesh& mesh
        );

        //- Constructor for gnemdFOAM
        pdBoundaries
        (
            Time& t,
            const polyMesh& mesh,
            pdCloud& molCloud
        );

        ~pdBoundaries();

        void setInitialConfig();

        void controlBeforeMove();
        void controlBeforeCollisions();
        void controlAfterCollisions();

        void calculateProps();

        //- update time information
        void updateTimeInfo();

        void outputResults();

        void checkPatchBoundaryModels(const polyMesh& mesh);
        void checkCyclicBoundaryModels(const polyMesh& mesh);


    // Member Functions

        const label& nPatchBoundaryModels() const;
        const label& nCyclicBoundaryModels() const;
        const label& nGeneralBoundaryModels() const;

        inline List< autoPtr<pdPatchBoundary> >& patchBoundaryModels();
        inline const List< autoPtr<pdPatchBoundary> >& patchBoundaryModels() const;
        inline const List<word>& patchBoundaryNames() const;
        inline const List<label>& patchBoundaryIds() const;
        inline const List<label>& patchToModelIds() const;

        inline List< autoPtr<pdCyclicBoundary> >& cyclicBoundaryModels();
        inline const List< autoPtr<pdCyclicBoundary> >& cyclicBoundaryModels() const;
        inline const List<word>& cyclicBoundaryNames() const;
        inline const List<label>& cyclicBoundaryIds() const;
        inline const List<label>& cyclicBoundaryToModelIds() const;


        inline List< autoPtr<pdGeneralBoundary> >& generalBoundaryModels();
        inline const List< autoPtr<pdGeneralBoundary> >& generalBoundaryModels() const;
        inline const List<word>& generalBoundaryNames() const;
        inline const List<label>& generalBoundaryIds() const;

    // Ostream Operator

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "pdBoundariesI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
