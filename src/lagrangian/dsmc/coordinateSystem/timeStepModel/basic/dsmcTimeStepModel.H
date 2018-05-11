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
    dsmcTimeStepModel

Description

SourceFiles
    dsmcTimeStepModelI.H
    dsmcTimeStepModel.C

\*---------------------------------------------------------------------------*/

#ifndef dsmcTimeStepModel_H
#define dsmcTimeStepModel_H

#include "fvMesh.H"
#include "volFields.H"
#include "dsmcParcel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class dsmcTimeStepModel Declaration
\*---------------------------------------------------------------------------*/


class dsmcTimeStepModel
{

protected:
    
    // Protected data
        
        //- Const reference to the mesh
        const fvMesh& mesh_;
        
        //- Number of real particles represented by a DSMC parcel
        volScalarField nParticles_;

        
        
private:

    // Private Member Functions
    
        //- Disallow default bitwise copy construct
        dsmcTimeStepModel(const dsmcTimeStepModel&);

        //- Disallow default bitwise assignment
        void operator=(const dsmcTimeStepModel&);



public:

    friend class dsmcParcel;

    //- Runtime type information
    TypeName("dsmcTimeStepModel");
    
    
    // Declare run-time constructor selection table
    #ifndef SWIG
    declareRunTimeSelectionTable
    (
        autoPtr,
        dsmcTimeStepModel,
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
        dsmcTimeStepModel
        (
            Time& t,
            const polyMesh& mesh,
            dsmcCloud& cloud
        );


    // Selectors	  
  
        //- Return a reference to the selected time-step model
        static autoPtr<dsmcTimeStepModel> New
        (
            Time& t,
            const polyMesh& mesh,
            dsmcCloud& cloud
        );
        
        
    // Destructor

        virtual ~dsmcTimeStepModel();


    // Member Functions

      // Access
      
        //- Return const access to the time-step value as given by the 
        //  system/controlDict dictionary
        inline scalar deltaTValueOrg() const;
        
        //- Return const access to the time-step field
        virtual inline const scalarField& deltaTField() const;
        
        //- Return const access to the time-step for a given cell
        virtual inline scalar deltaTValue(label celli = 0) const;
        
        //- Return const access to the number of real particles represented by
        //  a DSMC parcel
        virtual inline const volScalarField& nParticles() const;
        
        //- Return const access to the number of real particles represented by
        //  a DSMC parcel for a given cell
        virtual inline scalar nParticles(const label celli) const;
        
        //- Return const access to the number of real particles represented by
        //  a parcel for a given patch/face couple
        virtual inline scalar nParticles
        (
            const label patchi,
            const label facei
        ) const;
        
        //- Return non-const access to the number of real particles represented
        //  by a DSMC parcel for a given cell
        virtual inline scalar nParticles(const label celli);
        
        //- Return non-const access to the number of real particles represented
        //  by a parcel for a given patch/face couple
        virtual inline scalar nParticles
        (
            const label patchi,
            const label facei
        );
        
        
      // Check
      
        //- Check the validity of the time-step model inputs present in
        //  constant/dsmcProperties and assign inputs to data members if
        //  correct.
        virtual void checkTimeStepModelInputs() = 0;
      
      
      // Edit
      
        //- Update function
        virtual void update() = 0;
        
        //- Initialise nParticles field
        void initialisenParticles(const scalar value);
        

      // Write
      
        //- Write time-step model info
        virtual void writeTimeStepModelInfo() const = 0;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "dsmcTimeStepModelI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
