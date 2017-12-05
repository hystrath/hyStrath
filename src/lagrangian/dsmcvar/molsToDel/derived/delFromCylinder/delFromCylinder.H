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
    delFromCylinder

Description

SourceFiles
    delFromCylinderI.H
    delFromCylinder.C
    delFromCylinderIO.C

\*---------------------------------------------------------------------------*/

#ifndef delFromCylinder_H
#define delFromCylinder_H

#include "molsToDeleteModel.H"
#include "vector.H"
#include "polyMesh.H"
#include "volFields.H"
#include "dsmcCloud.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class delFromCylinder Declaration
\*---------------------------------------------------------------------------*/

class delFromCylinder
:
	public molsToDeleteModel
{

private:

    // Private data

        //- subDictionary containing the properties
        dictionary propsDict_;

        vector startPoint_;
        vector endPoint_;
        vector unitVector_;

        scalar radius_;

        labelList typeIds_;



    // Private Member Functions

        void findMolsToDel();

public:

    //- Runtime type information
    TypeName("delFromCylinder");

    // Constructors

        //- Construct from components
        delFromCylinder
		(
            dsmcCloud& cloud,
            const dictionary& dict
		);

    // Destructor

 		~delFromCylinder();


    // Member Functions

 
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
