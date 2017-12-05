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
    dsmcFullyCatalyticDiffuseWallPatch

Description
    Wall interaction setting microscopic velocity to a random one
    drawn from a Maxwellian distribution corresponding to a specified
    temperature

    SINGLE VALUE: temperature and velocity required in boundariesDict


SourceFiles
    dsmcFullyCatalyticDiffuseWallPatch.C


\*---------------------------------------------------------------------------*/

#ifndef dsmcFullyCatalyticDiffuseWallPatch_H
#define dsmcFullyCatalyticDiffuseWallPatch_H

#include "dsmcDiffuseWallPatch.H"
#include "dsmcCatalyticWallPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
             Class dsmcFullyCatalyticDiffuseWallPatch Declaration
\*---------------------------------------------------------------------------*/

class dsmcFullyCatalyticDiffuseWallPatch
:
	  public dsmcDiffuseWallPatch,
	  public dsmcCatalyticWallPatch
{

private:

    // Private data
    
        dictionary propsDict_;
        

protected:

    // Protected data
    

        
    // Protected member functions
      
        void setProperties();


public:

    //- Runtime type information
    TypeName("dsmcFullyCatalyticDiffuseWallPatch");

    
    // Constructors

        //- Construct from components
        dsmcFullyCatalyticDiffuseWallPatch
        (
            Time& t,
            const polyMesh& mesh,
            dsmcCloud& cloud,
            const dictionary& dict
        );

    // Destructor

        ~dsmcFullyCatalyticDiffuseWallPatch();


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
        
        void writeField();

        //- update properties from a modified dictionary
        void updateProperties(const dictionary&);



    // Access





};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
