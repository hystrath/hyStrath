/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016-2020 hyStrath
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of hyStrath, a derivative work of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    pdSpecularWallPatch

Description

     Specular wall boundary patch: reflecting molecules back in by reversing their
        normal velocity vector to the boundary

SourceFiles
    pdSpecularWallPatch.C


\*---------------------------------------------------------------------------*/

#ifndef pdSpecularWallPatch_H
#define pdSpecularWallPatch_H

#include "pdPatchBoundary.H"
#include "vector.H"
#include "polyMesh.H"
#include "volFields.H"
#include "pdCloud.H"
#include "distribution.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class pdSpecularWallPatch Declaration
\*---------------------------------------------------------------------------*/

class pdSpecularWallPatch
:
	public pdPatchBoundary
{

private:

    // Private data
        dictionary propsDict_;


public:

    //- Runtime type information
    TypeName("pdSpecularWallPatch");

    // Constructors

        //- Construct from components
        pdSpecularWallPatch
		(
            Time& t,
            const polyMesh& mesh,
            pdCloud& cloud,
            const dictionary& dict
		);

    // Destructor

 		~pdSpecularWallPatch();


    // Member Functions

        void initialConfiguration();

        //- calculate any required properties
        void calculateProperties();

        //- control the mols
        void controlParticle(pdParcel& p, pdParcel::trackingData& td);

        //- output data
        void output
        (
            const fileName& fixedPathName,
            const fileName& timePath
        );

        //- update properties from a modified dictionary
        void updateProperties(const dictionary&);



        // Access





};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
