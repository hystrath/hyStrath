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
    dsmcCartesian

Description

SourceFiles
    dsmcCartesianI.H
    dsmcCartesian.C

\*---------------------------------------------------------------------------*/

#ifndef dsmcCartesian_H
#define dsmcCartesian_H

#include "dsmcCoordinateSystem.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class dsmcCartesian Declaration
\*---------------------------------------------------------------------------*/


class dsmcCartesian
:
    public dsmcCoordinateSystem
{
    
    // Private data


    // Private Member Functions



public:

    //- Runtime type information
    TypeName("dsmcCartesian");
    
    
    // Constructors

        //- Constructor
        dsmcCartesian
        (
            Time& t,
            const polyMesh& mesh,
            dsmcCloud& cloud
        );


    // Destructor

        ~dsmcCartesian();   


    // Member Functions

      // Access
        
        
      // Check
      
        //- Check the validity of the coordinate system inputs present in
        //  constant/dsmcProperties and assign inputs to data members if
        //  correct. Pass init as true if used in the constructor of 
        //  the dsmcInitialise application
        void checkCoordinateSystemInputs(const bool init = false);
      
      
      // Edit
      
        //- Evolve function
        virtual void evolve();


      // Write
      
        //- Write coordinate system info
        virtual void writeCoordinateSystemInfo() const;
      

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "dsmcCartesianI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
