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

\*---------------------------------------------------------------------------*/

#ifndef gravitationalAccelerationControllerRadial_H
#define gravitationalAccelerationControllerRadial_H

#include "dsmcStateController.H"
#include "vector.H"
#include "fvMesh.H"
#include "volFields.H"
#include "dsmcCloud.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
Class gravitationalAccelerationControllerRadial Declaration
\*---------------------------------------------------------------------------*/

class gravitationalAccelerationControllerRadial
    :
        public dsmcStateController
{

    private:
        //- subDictionary containing the properties
        dictionary propsDict_;

        
        // use right hand grip rule to define rotational direction of acceleration
        // important to define well otherwise flow can be reversed
        vector rotationalAxis_;

        // midpoint between rotational start and end points
        vector m_;
        
        // define two reference vectors orthogonal to each other and the rotational axis
        vector rVx_;
        vector rVy_;
        
        // radius 
//         scalar radius_;
        
        // angle with respect to the reference vectors (important to define clockwise)
        scalar thetaStart_;
        scalar thetaEnd_;        
        
        scalar acc_;

        // Private Member Functions
        void setProperties();
        
        
        
    public:

        //- Runtime type information
        TypeName("gravitationalAccelerationControllerRadial");

        // Constructors
        //- Construct from components
        gravitationalAccelerationControllerRadial
            (
                Time& t,
                dsmcCloud& cloud,
                const dictionary& dict
            );

        // Destructor
        ~gravitationalAccelerationControllerRadial();

        // Member Functions
        void initialConfiguration();
        void controlParcelsBeforeMove();
        void controlParcelsBeforeCollisions();
        void controlParcelsAfterCollisions();

        //- calculate any required properties
        void calculateProperties();

        //- output data
        void output
            (
                const fileName& fixedPathName,
                const fileName& timePath
            );

        //- update properties from a modified dictionary
        void updateProperties(const dictionary&);
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
