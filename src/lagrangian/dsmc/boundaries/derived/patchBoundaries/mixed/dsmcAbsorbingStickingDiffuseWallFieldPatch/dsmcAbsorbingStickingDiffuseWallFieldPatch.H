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
    dsmcAbsorbingStickingDiffuseWallFieldPatch
    
Description    
    Wall interaction where a particle can be adsorbed according to a 
    probability. Othwerwise, its properties are drawn from a Maxwellian 
    distribution corresponding to a specified temperature
    
    FIELD VERSION: require boundaryU and boundaryT fields to be defined in
    your time directory.
    
SourceFiles
    dsmcAbsorbingStickingDiffuseWallFieldPatch.C
\*---------------------------------------------------------------------------*/

#ifndef dsmcAbsorbingStickingDiffuseWallFieldPatch_H
#define dsmcAbsorbingStickingDiffuseWallFieldPatch_H

#include "dsmcFieldPatchBoundary.H"
#include "dsmcAbsorbingWallPatch.H"
#include "dsmcStickingWallPatch.H"
#include "dsmcDiffuseWallPatch.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
          Class dsmcAbsorbingStickingDiffuseWallFieldPatch Declaration
\*---------------------------------------------------------------------------*/

class dsmcAbsorbingStickingDiffuseWallFieldPatch
:
    public dsmcFieldPatchBoundary,
    public dsmcDiffuseWallPatch,
    public dsmcAbsorbingWallPatch,
    public dsmcStickingWallPatch
{

private:

    // Private data
    
        
        
protected:        

    // Protected member functions
            
        

public:

    //- Runtime type information
    TypeName("dsmcAbsorbingStickingDiffuseWallFieldPatch");

    // Constructors

        //- Construct from components
        dsmcAbsorbingStickingDiffuseWallFieldPatch
        (
            Time& t,
            const polyMesh& mesh,
            dsmcCloud& cloud,
            const dictionary& dict
        );

    // Destructor

        ~dsmcAbsorbingStickingDiffuseWallFieldPatch();


    // Member Functions

        virtual void initialConfiguration();

        //- calculate any required properties
        virtual void calculateProperties();

        //- control the mols
        virtual void controlParticle
        (
            dsmcParcel& p,
            dsmcParcel::trackingData& td
        );

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
    
    
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
