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
    dsmcVariableTimeStepModel

Description

SourceFiles
    dsmcVariableTimeStepModelI.H
    dsmcVariableTimeStepModel.C

\*---------------------------------------------------------------------------*/

#ifndef dsmcVariableTimeStepModel_H
#define dsmcVariableTimeStepModel_H

#include "dsmcTimeStepModel.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class dsmcVariableTimeStepModel Declaration
\*---------------------------------------------------------------------------*/


class dsmcVariableTimeStepModel
:
    public dsmcTimeStepModel
{
    
    // Private data
    
      //- Const reference to the DSMC cloud
      const dsmcCloud& cloud_;
        
      //- Label of the reference cell for the variable time-step method
      //  described in [Wu, 2004]
      label refCell_;
      
      //- Variable time-step
      volScalarField deltaT_;
        

    // Private Member Functions
    
      //- Find the reference cell (cell with the smallest volume in the domain)
      //  and edit refCell_
      void findRefCell();
      
      //- Calculate number of equivalent particles field according to the
      //  variable time step method described in [Wu, 2004]
      //  Uses the cell volumes as an input
      void updatenParticles();
      
      //- Calculate the variable time step field according to the 
      //  variable time-step method described in [Wu, 2004]
      //  Uses the number of equivalent particles and the reference cell as an
      //  input
      void updateTimeStep();
      
      //- Update time step and nParticle fields
      //  when the mesh is changed
      void updateVariableTimeStepMethod();



public:

    //- Runtime type information
    TypeName("dsmcVariableTimeStepModel");
    
    // Constructors

        //- Constructor
        dsmcVariableTimeStepModel
        (
            Time& t,
            const polyMesh& mesh,
            dsmcCloud& cloud
        );


    // Destructor

        ~dsmcVariableTimeStepModel();


    // Member Functions

      // Access
      
        //- Return const access to the time-step field
        virtual inline const scalarField& deltaTField() const;
        
        //- Return const access to the time-step for a given cell
        virtual inline scalar deltaTValue(const label celli = 0) const;
        
        
      // Check
      
        //- Check the validity of the time-step model inputs present in
        //  constant/dsmcProperties and assign inputs to data members if
        //  correct.
        virtual void checkTimeStepModelInputs();
      
      
      // Edit
      
        //- Update function
        virtual void update();
        

      // Write
      
        //- Write time-step model info
        virtual void writeTimeStepModelInfo() const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "dsmcVariableTimeStepModelI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
