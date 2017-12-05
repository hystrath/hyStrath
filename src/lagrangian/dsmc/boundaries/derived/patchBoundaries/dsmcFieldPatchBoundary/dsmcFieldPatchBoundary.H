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
    dsmcFieldPatchBoundary

Description
    Wall interaction setting microscopic velocity to a random one
    drawn from a Maxwellian distribution corresponding to a specified
    temperature

    FIELD VERSION: require boundaryU and boundaryT fields to be defined in
    your time directory.

SourceFiles
    dsmcFieldPatchBoundary.C


\*---------------------------------------------------------------------------*/

#ifndef dsmcFieldPatchBoundary_H
#define dsmcFieldPatchBoundary_H

#include "dsmcPatchBoundary.H"

#include "dsmcCloud.H"
#include "vectorList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class dsmcFieldPatchBoundary Declaration
\*---------------------------------------------------------------------------*/

class dsmcFieldPatchBoundary
:
	  public virtual dsmcPatchBoundary
{

private:

    // Private data
    

protected:

    // Protected data        

        //- boundary temperature
        //volScalarField boundaryT_;
        scalarList boundaryT_;

        //- boundary velocity
        //volVectorField boundaryU_;
        vectorList boundaryU_;
        
        
    // Protected member functions
    
        //- Read the temperature and velocity fields on the patch
        void readPatchFields();
    
        

public:

    //- Runtime type information
    TypeName("dsmcFieldPatchBoundary");

    // Constructors

        //- Construct from components
        dsmcFieldPatchBoundary
	(
            Time& t,
            const polyMesh& mesh,
            dsmcCloud& cloud,
            const dictionary& dict
	);

    // Destructor

 		~dsmcFieldPatchBoundary();


    // Member Functions

        void initialConfiguration();

        //- calculate any required properties
        void calculateProperties();

        //- control the mols
        virtual void controlParticle
        (
            dsmcParcel& p,
            dsmcParcel::trackingData& td
        ) = 0;

        //- output data
        void output
        (
            const fileName& fixedPathName,
            const fileName& timePath
        );

        //- update properties from a modified dictionary
        void updateProperties(const dictionary&);


    // Access
    
      //- Return the patch temperature at the parcel
      scalar patchLocalTemperature(const dsmcParcel& p) const;
        
      //- Return the patch velocity at the parcel  
      const vector& patchLocalVelocity(const dsmcParcel& p) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
