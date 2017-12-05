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
    dsmcDiffuseWallPatch

Description
    Wall interaction setting microscopic velocity to a random one
    drawn from a Maxwellian distribution corresponding to a specified
    temperature

    SINGLE VALUE: temperature and velocity required in boundariesDict


SourceFiles
    dsmcDiffuseWallPatch.C


\*---------------------------------------------------------------------------*/

#ifndef dsmcDiffuseWallPatch_H
#define dsmcDiffuseWallPatch_H

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
                    Class dsmcDiffuseWallPatch Declaration
\*---------------------------------------------------------------------------*/

class dsmcDiffuseWallPatch
:
	  public virtual dsmcPatchBoundary
{

private:

    // Private data
    
        dictionary propsDict_;
        
        
protected: 

    // Protected data       
        
        //- depthAxis indicates the direction of depth for the medium
        //  default: y-axis (label = 1)
        label depthAxis_;
        
        //- Maximum coordinate (that of the ground level) along the depthAxis
        //  direction
        scalar maxDepth_;
        
        //- Length of the patch along the depthAxis taken as the difference
        //  between the max and min of the patch bounding box
        scalar lengthPatch_;
        
        //- formationLevelTemperature_ is the temperature to be set on the 
        //  formation surface of the porous medium. Linear temperature gradient
        //  between temperature_ that is the ground level temperature and 
        //  formationLevelTemperature_
        //  default: temperature_, i.e., no temperature gradient
        scalar formationLevelTemperature_;
        

    // Protected member functions
    
        virtual void setProperties();

        //- Perform a diffuse particle-wall reflection
        //  Local info on velocity and temperature are optional and overwrites
        //  default temperature_ and velocity_ attributs if provided
        void performDiffuseReflection
        (
            dsmcParcel& p,
            const scalar& localTemperature = 0,
            const vector& localVelocity = vector::zero
        );

        //- Return the temperature as a function of depth
        scalar getLocalTemperature(const scalar& depthPosition) const;
        

public:

    //- Runtime type information
    TypeName("dsmcDiffuseWallPatch");

    // Constructors

        //- Construct from components
        dsmcDiffuseWallPatch
        (
            Time& t,
            const polyMesh& mesh,
            dsmcCloud& cloud,
            const dictionary& dict
        );

    // Destructor

        ~dsmcDiffuseWallPatch();


    // Member Functions

        virtual void initialConfiguration();

        //- calculate any required properties
        virtual void calculateProperties();

        //- control the mols
        virtual void controlParticle(dsmcParcel& p, dsmcParcel::trackingData& td);
        
        //- output data
        virtual void output
        (
            const fileName& fixedPathName,
            const fileName& timePath
        );
        
        void writeField();

        //- update properties from a modified dictionary
        virtual void updateProperties(const dictionary&);


    // Access

        //- return const access to the depthAxis label
        const label& depthAxis() const
        {
            return depthAxis_;
        }
        
        //- return const access to the formation level temperature
        const scalar& formationLevelTemperature() const
        {
            return formationLevelTemperature_;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
