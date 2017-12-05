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
    dsmcSphericalPatchBoundary

Description
    Wall interaction setting microscopic velocity to a random one
    drawn from a Maxwellian distribution corresponding to a specified
    temperature

    FIELD VERSION: require boundaryU and boundaryT fields to be defined in
    your time directory.

SourceFiles
    dsmcSphericalPatchBoundary.C


\*---------------------------------------------------------------------------*/

#ifndef dsmcSphericalPatchBoundary_H
#define dsmcSphericalPatchBoundary_H

#include "dsmcPatchBoundary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class dsmcSphericalPatchBoundary Declaration
\*---------------------------------------------------------------------------*/

class dsmcSphericalPatchBoundary
:
	  public virtual dsmcPatchBoundary
{

private:

    // Private data
    
        dictionary propsDict_;
        

protected:

    // Protected data        

        //- Centre point, default: (0, 0, 0)
        const vector centrePoint_;
        
        //- Angular velocities, default: 0.0
        const scalar angularVelocityXPlane_;
        
        const scalar angularVelocityYPlane_;
        
        const scalar angularVelocityZPlane_;
        
        
    // Protected member functions
    
        vector wallVelocity(const dsmcParcel& p);
        
        

public:

    //- Runtime type information
    TypeName("dsmcSphericalPatchBoundary");

    // Constructors

        //- Construct from components
        dsmcSphericalPatchBoundary
		    (
            Time& t,
            const polyMesh& mesh,
            dsmcCloud& cloud,
            const dictionary& dict
		    );

    // Destructor

 		~dsmcSphericalPatchBoundary();


    // Member Functions

        void initialConfiguration();

        //- calculate any required properties
        void calculateProperties();

        //- control the mols
        virtual void controlParticle(dsmcParcel& p, dsmcParcel::trackingData& td) = 0;

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
