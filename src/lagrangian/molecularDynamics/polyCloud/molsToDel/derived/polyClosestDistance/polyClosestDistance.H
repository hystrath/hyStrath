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
    polyClosestDistance

Description
    Model which deletes polyMolecules from the mesh. Method should be used to delete
    polyMolecules using their id number. If oneSpecie = no, then all polyMolecules from
    the mesh are deleted.


SourceFiles
    polyClosestDistanceI.H
    polyClosestDistance.C
    polyClosestDistanceIO.C

\*---------------------------------------------------------------------------*/

#ifndef polyClosestDistance_H
#define polyClosestDistance_H

#include "polyMolsToDeleteModel.H"
#include "vector.H"
#include "polyMesh.H"
#include "volFields.H"
#include "polyMoleculeCloud.H"
#include "boundBox.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class polyClosestDistance Declaration
\*---------------------------------------------------------------------------*/

class polyClosestDistance
:
	public polyMolsToDeleteModel
{

private:

    // Private data

        //- subDictionary containing the properties
        dictionary propsDict_;
        
        scalar distance_;
//         List<label> molIds_;

        // mol id of molecules to be deleted
        label molToDeleteId_;
        
        // mol id of reference molecules
        label refMolId_;
        
        
    // Private Member Functions

        void checkBoundBox
        (
            boundBox& b,
            const vector& startPoint,
            const vector& endPoint
        );

        void findMolsToDel();

public:

    //- Runtime type information
    TypeName("polyClosestDistance");

    // Constructors

        //- Construct from components
        polyClosestDistance
		(
            polyMoleculeCloud& molCloud,
            const dictionary& dict
		);

    // Destructor

 		~polyClosestDistance();


    // Member Functions

 
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
