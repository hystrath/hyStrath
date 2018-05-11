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
    dsmcAbsorbingWallPatch

Description
    Wall interaction setting microscopic velocity to a random one
    drawn from a Maxwellian distribution corresponding to a specified
    temperature

    SINGLE VALUE: temperature and velocity required in boundariesDict


SourceFiles
    dsmcAbsorbingWallPatch.C


\*---------------------------------------------------------------------------*/

#ifndef dsmcAbsorbingWallPatch_H
#define dsmcAbsorbingWallPatch_H

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
                    Class dsmcAbsorbingWallPatch Declaration
\*---------------------------------------------------------------------------*/

class dsmcAbsorbingWallPatch
:
	  public virtual dsmcPatchBoundary
{

private:

    // Private data
    
        dictionary propsDict_;
        
        
protected:

    // Protected data  
    
        //- typeIds to be adsorbed
        labelList typeIds_;
        
        //- probability of absorption for each species in typeIds_
        scalarList absorptionProbs_;      
        
        //- Saturation limit per face [#]
        scalarList saturationLimit_;
        


    // Protected member functions
    
        virtual void setProperties();
        
        //- Read the number of absorbed particles on the patch
        void readPatchField();
        
        //- Check if absorption can be still operated on this patch/face
        bool isNotSaturated(const label patchi, const label facei);
        
        //- Absorb the particle
        void absorbParticle
        (
            const label patchi,
            const label facei,
            dsmcParcel::trackingData& td
        );
        
        //- Test the particle for absorption
        void testParticleForAbsorption
        (
            dsmcParcel& p,
            dsmcParcel::trackingData& td
        );


public:

    //- Runtime type information
    TypeName("dsmcAbsorbingWallPatch");

    // Constructors

        //- Construct from components
        dsmcAbsorbingWallPatch
        (
            Time& t,
            const polyMesh& mesh,
            dsmcCloud& cloud,
            const dictionary& dict
        );

    // Destructor

        ~dsmcAbsorbingWallPatch();


    // Member Functions

        virtual void initialConfiguration();

        //- calculate any required properties
        virtual void calculateProperties();

        //- control the mols
        virtual void controlParticle
        (
            dsmcParcel& p, 
            dsmcParcel::trackingData& td
        ) = 0;

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
    
        //- Return const access to the patch saturation limit
        const scalarList& saturationLimit() const
        {
            return saturationLimit_;
        }
        
        //- Return const access to the patch saturation limit for a given face
        const scalar& saturationLimit(const label& facei) const
        {
            return saturationLimit_[facei];
        }
        
        //- Return const access to all absorption probabilities
        const List<scalar>& absorptionProbs() const
        {
            return absorptionProbs_;
        }
        
        //- Return const access to a molecule absorption probability
        //  Given by typeId
        scalar absorptionProb(const label& typeId) const
        {
            const label& iD = findIndex(typeIds_, typeId);
            
            if (iD != -1)
            {
                return absorptionProbs_[iD];
            }
            
            FatalErrorIn
            (
                "dsmcAbsorbingDiffuseWallPatch::const scalar&\
                     absorptionProb(const label& typeId) const"
            )
                << "typeId is not listed."
                << exit(FatalError);
            
            return 0;
        }

        
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
