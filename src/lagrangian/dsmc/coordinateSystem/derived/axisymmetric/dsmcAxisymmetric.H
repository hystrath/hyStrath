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
    dsmcAxisymmetric

Description

SourceFiles
    dsmcAxisymmetricI.H
    dsmcAxisymmetric.C

\*---------------------------------------------------------------------------*/

#ifndef dsmcAxisymmetric_H
#define dsmcAxisymmetric_H

#include "dsmcCoordinateSystem.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class dsmcAxisymmetric Declaration
\*---------------------------------------------------------------------------*/


class dsmcAxisymmetric
:
    public dsmcCoordinateSystem
{
    
    // Private data

        //- Reference to the DSMC cloud
        dsmcCloud& cloud_;
        
        //- Axis of revolution
        //  Default: x-axis, i.e., component 0
        label revolutionAxis_;
        
        //- Polar axis
        //  Default: (revolutionAxis_ + 1)%3
        label polarAxis_;
        
        //- Angular coordinate
        //  Default: (revolutionAxis_ + 2)%3
        label angularCoordinate_;
        
        //- Radial extent of the domain
        scalar radialExtent_;
        
        //- Maximum radial scaling factor
        scalar maxRWF_;
        
        //- Radial weighting method employed
        //  Default: "cell", other options: "particle", "mixed"
        word rWMethod_;
        
        //- Radial weighting factor field
        volScalarField RWF_;
        

    // Private Member Functions
    
        //- Clone or delete particles based on their old and new RWF
        void axisymmetricWeighting();
        
        //- Update the radial weighting factor field
        void updateRWF();
        
        //- Write axisymmetric info
        void writeAxisymmetricInfo() const;
        
        //- Return the radial extent of the domain
        inline scalar radialExtent() const;
        
        //- Return the maximum value for RWF
        inline scalar maxRWF() const;



public:

    //- Runtime type information
    TypeName("dsmcAxisymmetric");
    
    // Constructors

        //- Constructor
        dsmcAxisymmetric
        (
            Time& t,
            const polyMesh& mesh,
            dsmcCloud& cloud
        );


    // Destructor

        ~dsmcAxisymmetric();


    // Member Functions

      // Access
      
        //- Return const access to the number of real particles represented by
        //  a DSMC parcel
        virtual inline const volScalarField& nParticles() const;
        
        //- Return const access to the number of real particles represented by
        //  a DSMC parcel for a given cell
        //  Locally in the code, the radialWeightingMethod (RWM) can be
        //  changed from "cell"- to "particle"-based if rWMethod is set to
        //  mixed and the method call uses mixedRWM = true (e.g., nTC)
        virtual inline scalar nParticles
        (
            const label celli,
            const bool mixedRWMethod = false
        ) const;
        
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
        
        //- Return the local RWF for a given cell
        //  Locally in the code, the radialWeightingMethod (RWM) can be
        //  changed from "cell"- to "particle"-based if rWMethod is set to
        //  mixed and the method call uses mixedRWM = true (e.g., nTC)
        virtual inline scalar RWF
        (
            const label celli,
            const bool mixedRWMethod = false
        ) const;
        
        //- Return the local RWF for a given face
        virtual inline scalar pRWF
        (
            const label patchi,
            const label facei
        ) const;
        
        
      // Check
      
        //- Check the validity of the coordinate system inputs present in
        //  constant/dsmcProperties and assign inputs to data members if
        //  correct. Pass init as true if used in the constructor of 
        //  the dsmcInitialise application
        virtual void checkCoordinateSystemInputs(const bool init = false);
      
      
      // Edit
      
        //- Evolve function
        virtual void evolve();
        
        //- Recalculate and return the radial weighting factor for a given 
        //  patch/face
        virtual scalar recalculatepRWF
        (
            const label patchI,
            const label faceI
        ) const;
        
        //- Recalculate and return the radial weighting factor for a cell
        //  Locally in the code, the radialWeightingMethod (rWMethod) can be
        //  changed from "cell"- to "particle"-based if rWMethod is set to
        //  mixed and the method call uses mixedRWM = true (e.g., nTC)
        virtual scalar recalculateRWF
        (
            const label cellI, 
            const bool mixedRWMethod = false
        ) const;


      // Write
      
        //- Write coordinate system info
        virtual void writeCoordinateSystemInfo() const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "dsmcAxisymmetricI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
