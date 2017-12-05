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
    adaptiveLangevin

Description
    Written by Matthew Borg and Xiaocheng Shang


SourceFiles
    adaptiveLangevinI.H
    adaptiveLangevin.C
    adaptiveLangevinIO.C

\*---------------------------------------------------------------------------*/

#ifndef adaptiveLangevin_H
#define adaptiveLangevin_H

#include "polyIntegrator.H"
#include "vector.H"
#include "polyMesh.H"
#include "volFields.H"
#include "polyMoleculeCloud.H"
// #include "boundedBox.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class adaptiveLangevin Declaration
\*---------------------------------------------------------------------------*/

class adaptiveLangevin
:
    public polyIntegrator
{

private:
    
        //- subDictionary containing the properties
        dictionary propsDict_;
    
        scalar deltaT_;
        
        scalar KE1_;
        
        scalar xi_;
        
        scalar invMu_;
        
        scalar gamma_;
        
        scalar T_;
        
        scalar Nd_;
        
        scalar kB_;
        
        scalar N_;
        
        scalar sigma_;
        
        DynamicList<scalar> kE1_;
        DynamicList<scalar> kE2_;
        DynamicList<scalar> angularKE1_;
        DynamicList<scalar> angularKE2_;
        
        
    // Private data
        
        void calculateKE(scalar& kE, scalar& angularKeSum);
        
        void calculateXi();
       
        void updateHalfVelocity();
        
        void write();       
        
        
public:

    //- Runtime type information
    TypeName("adaptiveLangevin");

    // Constructors

        //- Construct from components
        adaptiveLangevin
        (
            Time& t,
            polyMoleculeCloud& molCloud,
            const dictionary& dict
        );

    // Destructor

        ~adaptiveLangevin();


    // Member Functions

        void init();
        
        void evolve();
        


};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
