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
    polyCyclicBoundary

Description

    Used for applying instantaneous actions on molecules colliding with a cyclic 
    boundary.

SourceFiles
    polyCyclicBoundaryI.H
    polyCyclicBoundary.C
    polyCyclicBoundaryIO.C

\*---------------------------------------------------------------------------*/

#ifndef polyCyclicBoundary_H
#define polyCyclicBoundary_H

#include "IOdictionary.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "vector.H"
#include "volFields.H"
#include "polyMesh.H"
#include "polyMolecule.H"
#include "Time.H"
#include "writeTimeData.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class polyCyclicBoundary Declaration
\*---------------------------------------------------------------------------*/

class polyCyclicBoundary
{

protected:

    // Protected data

        const fvMesh& mesh_;

        polyMoleculeCloud& molCloud_;

        Time& time_;

        //- name of this patch
        word patchName_;
        
        //- id on boundaryMesh
        label patchId_;
        
        //- writing out properties
        bool writeInTimeDir_;
        bool writeInCase_;

    // Private Member Functions

public:

    //- Runtime type information
    TypeName("polyCyclicBoundary");

    // Declare runtime constructor selection table
        declareRunTimeSelectionTable
        (
            autoPtr,
            polyCyclicBoundary,
            dictionary,
            (
                Time& t,
                const polyMesh& mesh,
                polyMoleculeCloud& molCloud,
                const dictionary& dict
            ),
            (t, mesh, molCloud, dict)
        );
   
    // Constructors

        //- Construct from components
        polyCyclicBoundary
        (
            Time& t,
            const polyMesh& mesh,
            polyMoleculeCloud& molCloud,
            const dictionary& dict
        );


    // Selectors

        static autoPtr<polyCyclicBoundary> New
        (
            Time& t,
            const polyMesh& mesh,
            polyMoleculeCloud& molCloud,
            const dictionary& dict
        );


    // Destructor

        virtual ~polyCyclicBoundary();


    // Member Functions

        //- update properties from a modified dictionary
        virtual void updateProperties(const dictionary&) = 0;

        void updateBoundaryProperties(const dictionary&);

        virtual void calculateProperties() = 0;

        //- create an initial configuration
        virtual void initialConfiguration() = 0;

        //- control the molecule
        virtual void controlMol
        (
            polyMolecule& mol,
            polyMolecule::trackingData& td,
            const label& patchi
        ) = 0;

        virtual void controlAfterMove() = 0;

        //- output data
        virtual void output
        (
            const fileName& fixedPathName,
            const fileName& timePath
        ) = 0;
        
	// Access

		const word& patchName() const;
		const label& patchId() const;
		const bool& writeInTimeDir() const;
		const bool& writeInCase() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
