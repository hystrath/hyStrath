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
    porousMeasurements

Description

SourceFiles
    porousMeasurementsI.H
    porousMeasurements.C

\*---------------------------------------------------------------------------*/

#ifndef porousMeasurements_H
#define porousMeasurements_H

#include "fvMesh.H"
#include "dsmcParcel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class porousMeasurements Declaration
\*---------------------------------------------------------------------------*/


class porousMeasurements
{

protected:
    
    // Protected data
    
        
        
private:

    // Private Member Functions
    
        //- Disallow default bitwise copy construct
        porousMeasurements(const porousMeasurements&);

        //- Disallow default bitwise assignment
        void operator=(const porousMeasurements&);



public:

    friend class dsmcParcel;

    //- Runtime type information
    TypeName("porousMeasurements");
    
    
    // Declare run-time constructor selection table
    #ifndef SWIG
    declareRunTimeSelectionTable
    (
        autoPtr,
        porousMeasurements,
        fvMesh,
        (
            Time& t,
            const polyMesh& mesh,
            dsmcCloud& cloud
        ),
        (t, mesh, cloud)
    );
    #endif
    
    
    // Constructors

        //- Constructor
        porousMeasurements
        (
            Time& t,
            const polyMesh& mesh,
            dsmcCloud& cloud
        );
        

    // Selectors	  
  
        //- Return a reference to the selected porous measurements model
        static autoPtr<porousMeasurements> New
        (
            Time& t,
            const polyMesh& mesh,
            dsmcCloud& cloud
        );
        
        
    // Destructor

        virtual ~porousMeasurements();


    // Member Functions

      // Access
      
        
        
      // Check
      
        //- Check the validity of the porous measurement inputs present in
        //  constant/dsmcProperties and assign inputs to data members if
        //  correct.
        virtual void checkPorousMeasurementsInputs() = 0;
      
      
      // Edit
      
        //- Update function
        virtual void update() = 0;
        
        //- Effect of a specular interaction on porous medium measurements
        virtual void specularInteraction(dsmcParcel& p, const vector& nw) = 0;
        
        //- Effect of a diffuse interaction on porous medium measurements
        virtual void diffuseInteraction(dsmcParcel& p) = 0;
        
        //- Effect of a cyclic membrane interaction on porous medium 
        //  measurements
        virtual void cyclicMembraneInteraction
        (
            dsmcParcel& p,
            const vector& orgPos
        ) {}
        
        //- Effect of a deletion patch on porous medium measurements
        virtual void deletionInteraction(dsmcParcel& p, const label patchId) {}
        
        //- Tracer selection for porous medium measurements
        virtual void additionInteraction(dsmcParcel& p, const label patchId) {}


      // Write
      
        //- Write porous measurements info
        virtual void writePorousMeasurementsInfo() const = 0;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "porousMeasurementsI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
