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
    dsmcDiffuseSpecularWallPatch

Description



SourceFiles
    dsmcDiffuseSpecularWallPatch.C


\*---------------------------------------------------------------------------*/

#ifndef dsmcDiffuseSpecularWallPatch_H
#define dsmcDiffuseSpecularWallPatch_H

#include "dsmcSpecularWallPatch.H"
#include "dsmcDiffuseWallPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
               Class dsmcDiffuseSpecularWallPatch Declaration
\*---------------------------------------------------------------------------*/

class dsmcDiffuseSpecularWallPatch
:
	  public dsmcSpecularWallPatch,
	  public dsmcDiffuseWallPatch
{

private:

    // Private data
        
        dictionary propsDict_;
        

protected:

    // Protected data
    
        //- Fraction of wall interactions that are diffuse
        scalar diffuseFraction_;
    
    
    // Protected member functions
    
        void setProperties();

public:

    //- Runtime type information
    TypeName("dsmcDiffuseSpecularWallPatch");

    // Constructors

        //- Construct from components
        dsmcDiffuseSpecularWallPatch
	      (
            Time& t,
            const polyMesh& mesh,
            dsmcCloud& cloud,
            const dictionary& dict
	      );

    // Destructor

      ~dsmcDiffuseSpecularWallPatch();


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

        //- Return const access to the diffuse fraction
        const scalar& diffuseFraction() const
        {
            return diffuseFraction_;
        }

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
