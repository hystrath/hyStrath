/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2009-2010 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.
Class
    Foam::LarsenBorgnakkeVariableHardSphereNoTv

Description
    Variable Hard Sphere BinaryCollision Model with Larsen Borgnakke internal
    energy redistribution.  Based on the INELRS subroutine in Bird's DSMC0R.FOR

\*---------------------------------------------------------------------------*/

#ifndef LarsenBorgnakkeVariableHardSphereNoTv_H
#define LarsenBorgnakkeVariableHardSphereNoTv_H

#include "BinaryCollisionModel.H"
#include "dsmcCloud.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                      Class LarsenBorgnakkeVariableHardSphereNoTv Declaration
\*---------------------------------------------------------------------------*/

//template<class CloudType>
class LarsenBorgnakkeVariableHardSphereNoTv
:
    public BinaryCollisionModel
{
    // Private data
    
    //- The coefficients dictionary
        const dictionary coeffDict_;

        //- Reference temperature
        const scalar Tref_;

        //- Rotational relaxation collision number
        const scalar rotationalRelaxationCollisionNumber_;

    // Private Member Functions


public:

    //- Runtime type information
    TypeName("LarsenBorgnakkeVariableHardSphereNoTv");


    // Constructors

        //- Construct from dictionary
        LarsenBorgnakkeVariableHardSphereNoTv
        (
            const dictionary& dict,
            dsmcCloud& cloud
        );


    // Destructor
    ~LarsenBorgnakkeVariableHardSphereNoTv();


    // Member Functions
    
    //- Return the coefficients dictionary
        const dictionary& coeffDict() const;

        //- Flag to indicate whether model activates collision model
        bool active() const;

        //- Return the collision cross section * relative velocity product
        scalar sigmaTcR
        (
            const dsmcParcel& pP,
            const dsmcParcel& pQ
        ) const;

        //- Apply collision
        void collide
        (
            dsmcParcel& pP,
            dsmcParcel& pQ
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// #ifdef NoRepository
// #   include "LarsenBorgnakkeVariableHardSphereNoTv.C"
// #endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
