/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2005 OpenCFD Ltd.
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
    dsmcGlobalPorousMediumMeasurements

Description

SourceFiles
    dsmcGlobalPorousMediumMeasurementsI.H
    dsmcGlobalPorousMediumMeasurements.C

\*---------------------------------------------------------------------------*/

#ifndef dsmcGlobalPorousMediumMeasurements_H
#define dsmcGlobalPorousMediumMeasurements_H

#include "porousMeasurements.H"
#include "volFields.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
            Class dsmcGlobalPorousMediumMeasurements Declaration
\*---------------------------------------------------------------------------*/


class dsmcGlobalPorousMediumMeasurements
:
    public porousMeasurements
{
    
    // Private data

        //- Const reference to the mesh
        const fvMesh& mesh_;
        
        //- Reference to the DSMC cloud
        dsmcCloud& cloud_;
        
        //- Indicate the fraction of the dsmcCloud that should be used
        scalar trackingProbability_;
        
        //- Map associating inflow patches with relevant outflow patches
        wordHashSet tracerPatchNamesMap_;
        
        //- Dimensionality of the problem (i.e., 1: 1D, 2: 2D or 3: 3D)
        label dimensionality_;
        
        //- Average flow transit time for all tracked parcels recorded 
        scalar mediumTransitTime_;
        
        //- Average total distance travelled by all tracked parcels recorded
        scalar mediumTotalDistanceTravelled_;
        
        //- Average tortuosity measured from all tracked parcels recorded
        scalar mediumTortuosity_;
        
        //- Number of tracked parcels exiting through wrong patches
        label nLooping_;
        
        //- Time spent by the tracked parcels referred to as 'looping'
        scalar timeLooping_;
        

    // Private Member Functions
    
        //- Update the medium properties on a deletion patch
        void updateMediumPropertiesMeasurement
        (
            dsmcParcel& p,
            const label& delPatchId
        );
        
        //- Update the medium properties on a cyclic patch
        void updateMediumPropertiesMeasurement_cyclic
        (
            dsmcParcel& p,
            const label& neiPatchId,
            const label& orgPatchId,
            const vector& orgPosition
        );
        
        //- Write instantaneous porous measurements info
        void writeGlobalPorousMeasurementsInfo() const;



public:

    //- Runtime type information
    TypeName("dsmcGlobalPorousMediumMeasurements");
    
    // Constructors

        //- Constructor
        dsmcGlobalPorousMediumMeasurements
        (
            Time& t,
            const polyMesh& mesh,
            dsmcCloud& cloud
        );


    // Destructor

        ~dsmcGlobalPorousMediumMeasurements();


    // Member Functions

      // Access
      
        
        
      // Check
      
        //- Check the validity of the porous measurement inputs present in
        //  constant/dsmcProperties and assign inputs to data members if
        //  correct.
        virtual void checkPorousMeasurementsInputs();
      
      
      // Edit
      
        //- Update function
        virtual void update();
        
        //- Effect of a specular interaction on porous medium measurements
        virtual void specularInteraction(dsmcParcel& p, const vector& nw);
        
        //- Effect of a diffuse interaction on porous medium measurements
        virtual void diffuseInteraction(dsmcParcel& p);
        
        //- Effect of a deletion patch on porous medium measurements
        virtual void deletionInteraction(dsmcParcel& p, const label patchId);
        
        //- Tracer selection for porous medium measurements
        virtual void additionInteraction(dsmcParcel& p, const label patchId);
        

      // Write
      
        //- Write porous measurements info
        virtual void writePorousMeasurementsInfo() const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "dsmcGlobalPorousMediumMeasurementsI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
