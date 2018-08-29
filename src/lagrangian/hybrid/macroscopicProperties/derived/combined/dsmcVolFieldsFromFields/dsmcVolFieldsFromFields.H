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
    dsmcVolFieldsFromFields

Description
    Geometric field for temperature, density, velocity, pressure.

SourceFiles
    dsmcVolFieldsFromFields.C

\*---------------------------------------------------------------------------*/

#ifndef dsmcVolFieldsFromFields_H
#define dsmcVolFieldsFromFields_H

#include "dsmcVolFields.H"
#include "vector.H"
#include "polyMesh.H"
#include "volFields.H"
#include "dsmcCloud.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class dsmcVolFieldsFromFields Declaration
\*---------------------------------------------------------------------------*/

class dsmcVolFieldsFromFields
:
	public dsmcVolFields
{

protected:

    // Protected data
    
      volVectorField heatFluxTraVector_;
      volVectorField heatFluxRotVector_;
      volVectorField heatFluxVibVector_;
      
      scalarField eRotU_;
      scalarField eRotV_;
      scalarField eRotW_;
      scalarField eVibU_;
      scalarField eVibV_;
      scalarField eVibW_;
      scalarField eVib_;


public:

    //- Runtime type information
    TypeName("dsmcVolFieldsFromFields");

    // Constructors

      //- Construct from components
      dsmcVolFieldsFromFields
      (
          Time& t,
          const polyMesh& mesh,
          dsmcCloud& cloud,
          const dictionary& dict
      );

    // Destructor

    ~dsmcVolFieldsFromFields();


    // Member Functions

      //- calculate field
      void calculateField();
      
      //- reset fields if mesh has changed
      void resetField();

      //- write field
      void writeField();

      void updateProperties(const dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

