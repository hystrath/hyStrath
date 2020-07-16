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
    dsmcAdiabaticWallRotationPatch

Description
    Wall interaction setting microscopic velocity to a random one
    drawn from a Maxwellian distribution corresponding to a specified
    temperature for a specified fraction of collisions, and reversing
    the wall-normal component of the particle velocity for the
    remainder.

    FIELD VERSION: require boundaryU and boundaryT fields to be defined in
    your time directory.

SourceFiles
    dsmcAdiabaticWallRotationPatch.C


\*---------------------------------------------------------------------------*/

#ifndef dsmcAdiabaticWallRotationPatch_H
#define dsmcAdiabaticWallRotationPatch_H

#include "dsmcPatchBoundary.H"
#include "vector.H"
#include "polyMesh.H"
#include "volFields.H"
#include "dsmcCloud.H"
#include "distribution.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class dsmcAdiabaticWallRotationPatch Declaration
\*---------------------------------------------------------------------------*/

class dsmcAdiabaticWallRotationPatch
:
	public dsmcPatchBoundary
{

private:

    // Private data
        dictionary propsDict_;

        vector centrePoint_;

        scalar wallVelocity_;

        scalar referenceTemperature_;

        vector rotationAxis_;

public:

    //- Runtime type information
    TypeName("dsmcAdiabaticWallRotationPatch");

    // Constructors

        //- Construct from components
        dsmcAdiabaticWallRotationPatch
        (
            Time& t,
            const polyMesh& mesh,
            dsmcCloud& cloud,
            const dictionary& dict
        );

    // Destructor

        ~dsmcAdiabaticWallRotationPatch();


    // Member Functions

        void initialConfiguration();

        //- calculate any required properties
        void calculateProperties();

        //- control the mols
        void controlParticle(dsmcParcel& p, dsmcParcel::trackingData& td);

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
