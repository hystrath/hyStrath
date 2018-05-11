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
    dsmcInstantaneousPorousMediumMeasurements

Description

SourceFiles
    dsmcInstantaneousPorousMediumMeasurementsI.H
    dsmcInstantaneousPorousMediumMeasurements.C

\*---------------------------------------------------------------------------*/

#ifndef dsmcInstantaneousPorousMediumMeasurements_H
#define dsmcInstantaneousPorousMediumMeasurements_H

#include "porousMeasurements.H"
#include "volFields.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
          Class dsmcInstantaneousPorousMediumMeasurements Declaration
\*---------------------------------------------------------------------------*/


class dsmcInstantaneousPorousMediumMeasurements
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
        
        //- Dimensionality of the problem (i.e., 1: 1D)
        label dimensionality_;
        
        //- Initial simulation time
        scalar initialSimulationTime_;
        
        //- Average total distance travelled by all tracked parcels recorded
        scalar mediumTotalDistanceTravelled_;
        

    // Private Member Functions
    
        //- Write instantaneous porous measurements info
        void writeInstantaneousPorousMeasurementsInfo() const;



public:

    //- Runtime type information
    TypeName("dsmcInstantaneousPorousMediumMeasurements");
    
    // Constructors

        //- Constructor
        dsmcInstantaneousPorousMediumMeasurements
        (
            Time& t,
            const polyMesh& mesh,
            dsmcCloud& cloud
        );


    // Destructor

        ~dsmcInstantaneousPorousMediumMeasurements();


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
        
        //- Effect of a cyclic membrane interaction on porous medium 
        //  measurements
        virtual void cyclicMembraneInteraction
        (
            dsmcParcel& p,
            const vector& orgPos
        );
        

      // Write
      
        //- Write porous measurements info
        virtual void writePorousMeasurementsInfo() const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "dsmcInstantaneousPorousMediumMeasurementsI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
