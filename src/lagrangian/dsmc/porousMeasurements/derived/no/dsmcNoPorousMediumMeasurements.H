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
    dsmcNoPorousMediumMeasurements

Description

SourceFiles
    dsmcNoPorousMediumMeasurementsI.H
    dsmcNoPorousMediumMeasurements.C

\*---------------------------------------------------------------------------*/

#ifndef dsmcNoPorousMediumMeasurements_H
#define dsmcNoPorousMediumMeasurements_H

#include "porousMeasurements.H"
#include "volFields.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
          Class dsmcNoPorousMediumMeasurements Declaration
\*---------------------------------------------------------------------------*/


class dsmcNoPorousMediumMeasurements
:
    public porousMeasurements
{
    
    // Private data
        

    // Private Member Functions



public:

    //- Runtime type information
    TypeName("dsmcNoPorousMediumMeasurements");
    
    // Constructors

        //- Constructor
        dsmcNoPorousMediumMeasurements
        (
            Time& t,
            const polyMesh& mesh,
            dsmcCloud& cloud
        );


    // Destructor

        ~dsmcNoPorousMediumMeasurements();


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
        

      // Write
      
        //- Write porous measurements info
        virtual void writePorousMeasurementsInfo() const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "dsmcNoPorousMediumMeasurementsI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
