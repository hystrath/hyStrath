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
    dsmcConstantTimeStepModel

Description

SourceFiles
    dsmcConstantTimeStepModelI.H
    dsmcConstantTimeStepModel.C

\*---------------------------------------------------------------------------*/

#ifndef dsmcConstantTimeStepModel_H
#define dsmcConstantTimeStepModel_H

#include "dsmcTimeStepModel.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class dsmcConstantTimeStepModel Declaration
\*---------------------------------------------------------------------------*/


class dsmcConstantTimeStepModel
:
    public dsmcTimeStepModel
{
    
    // Private data
    
    

    // Private Member Functions



public:

    //- Runtime type information
    TypeName("dsmcConstantTimeStepModel");
    
    // Constructors

        //- Constructor
        dsmcConstantTimeStepModel
        (
            Time& t,
            const polyMesh& mesh,
            dsmcCloud& cloud
        );


    // Destructor

        ~dsmcConstantTimeStepModel();


    // Member Functions

      // Access
      
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

#include "dsmcConstantTimeStepModelI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
