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
    dsmcMixedDiffuseSpecularSphericalPatch

Description
    Wall interaction setting microscopic velocity to a random one
    drawn from a Maxwellian distribution corresponding to a specified
    temperature for a specified fraction of collisions, and reversing
    the wall-normal component of the particle velocity for the
    remainder.

    FIELD VERSION: require boundaryU and boundaryT fields to be defined in
    your time directory.

SourceFiles
    dsmcMixedDiffuseSpecularSphericalPatch.C


\*---------------------------------------------------------------------------*/

#ifndef dsmcMixedDiffuseSpecularSphericalPatch_H
#define dsmcMixedDiffuseSpecularSphericalPatch_H

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
                           Class dsmcMixedDiffuseSpecularSphericalPatch Declaration
\*---------------------------------------------------------------------------*/

class dsmcMixedDiffuseSpecularSphericalPatch
:
	public dsmcPatchBoundary
{

private:

    // Private data
        dictionary propsDict_;

        //- boundary temperature
//         volScalarField boundaryT_;

        //- Fraction of wall interactions that are diffuse
        scalar diffuseFraction_;

        vector centrePoint_;

        scalar angularVelocityXPlane_;
        
        scalar angularVelocityYPlane_;
        
        scalar angularVelocityZPlane_;
        
        scalar wallTemperature_;
        
public:

    //- Runtime type information
    TypeName("dsmcMixedDiffuseSpecularSphericalPatch");

    // Constructors

        //- Construct from components
        dsmcMixedDiffuseSpecularSphericalPatch
        (
            Time& t,
            const polyMesh& mesh,
            dsmcCloud& cloud,
            const dictionary& dict
        );

    // Destructor

        ~dsmcMixedDiffuseSpecularSphericalPatch();


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
