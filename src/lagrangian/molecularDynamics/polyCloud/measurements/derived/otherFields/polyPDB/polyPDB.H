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
    polyPDB


Description
    to write out a pdb file for creating movies

SourceFiles
    polyPDBI.H
    polyPDB.C
    polyPDBIO.C

\*---------------------------------------------------------------------------*/

#ifndef polyPDB_H
#define polyPDB_H

#include "polyField.H"
#include "vector.H"
#include "polyMesh.H"
#include "volFields.H"
#include "polyMoleculeCloud.H"
#include "distribution.H"
#include "boundedBox.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class polyPDB Declaration
\*---------------------------------------------------------------------------*/

class polyPDB
:
	public polyField
{

private:

    // Private data

        //- subDictionary containing the properties
        dictionary propsDict_;

        propertyField fields_;

        List<label> molIds_;

        List<word> excludeSites_;

        word fieldName_;

        label n_;

        List<label> minLimit_;
        List<label> maxLimit_;

        label iteration_;

        
        word option_;
        
        // zone        
        
        word regionName_;
        label regionId_;

        // bound boxes
        
        List<boundedBox> boxes_;
        
        label timeIndex_;
        label nSteps_;
        
        bool variableMols_;
        label nSiteEstimate_;
        vector rDummy_;
        word molOption_;
        label nSitesMol_;
        
        scalar startTime_;
        scalar endTime_;
        scalar accumulatedTime_;
        scalar deltaT_;
    
        bool writeFirstTimeStep_;
        
    // Private Member Functions

            void write();

            void writeInZone(List<labelField>& molIds, List<vectorField>& sites);

            void writeInMesh(List<labelField>& molIds, List<vectorField>& sites);
            
            void writeInBoundBox(List<labelField>& molIds, List<vectorField>& sites);            

            
public:

    //- Runtime type information
    TypeName("polyPDB");

    // Constructors

        //- Construct from components
        polyPDB
		(
            Time& t,
            const polyMesh& mesh,
            polyMoleculeCloud& molCloud,
            const dictionary& dict
		);

    // Destructor

 		~polyPDB();

    // Member Functions

        //- create initial field
        void createField();

        //- calculate field
        void calculateField();

        void measureDuringForceComputation
        (
            polyMolecule* molI,
            polyMolecule* molJ
        );

        void measureDuringForceComputationSite
        (
            polyMolecule* molI,
            polyMolecule* molJ,
            label sI,
            label sJ
        );

        const propertyField& fields() const;

        //- write field
        void writeField();
}; 		


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
