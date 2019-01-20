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
    Foam::LarsenBorgnakkeVariableHardSphere

Description
    Variable Hard Sphere BinaryCollision Model with Larsen Borgnakke internal
    energy redistribution. Based on the INELRS subroutine in Bird's DSMC0R.FOR

\*---------------------------------------------------------------------------*/

#ifndef LarsenBorgnakkeVariableHardSphere_H
#define LarsenBorgnakkeVariableHardSphere_H

#include "VariableHardSphere.H"
#include "dsmcCloud.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
               Class LarsenBorgnakkeVariableHardSphere Declaration
\*---------------------------------------------------------------------------*/

class LarsenBorgnakkeVariableHardSphere
:
    public VariableHardSphere
{
    // Private data
    
        //- The coefficients dictionary
        const dictionary coeffDictLB_;

        //- Rotational relaxation collision number
        const scalar rotationalRelaxationCollisionNumber_;
        
        //- Vibrational relaxation collision number
        //  Default: 0, i.e., non constant (see equation 2, Bird 2010)
        //  If a value is given in the dsmcProperties dict, then it is affected
        //  to vibrationalRelaxationCollisionNumber_
        const scalar vibrationalRelaxationCollisionNumber_;
        
        //- Formulation for the inverse vibrational collision number
        //  that is: 
        //  pre-2008: 1/Zv using the collision temperature -> 0
        //  2008: 1/Zv using the overall temperature -> 1
        //  (Bird 2008 RGD "A Comparison of Collision-Based and Temp-Based...")
        //  post-2008: 1/(5*Zv) using the collision temperature -> 2
        label invZvFormulation_;
        
        //- Electronic relaxation collision number
        const scalar electronicRelaxationCollisionNumber_;


    // Private Member Functions


public:

    //- Runtime type information
    TypeName("LarsenBorgnakkeVariableHardSphere");


    // Constructors

        //- Construct from dictionary
        LarsenBorgnakkeVariableHardSphere
        (
            const dictionary& dict,
            dsmcCloud& cloud
        );


    // Destructor
    ~LarsenBorgnakkeVariableHardSphere();


    // Member Functions
    
        //- Return the vibrational relaxation number
        inline const scalar& vibrationalRelaxationNumber() const
        {
            return vibrationalRelaxationCollisionNumber_;
        }
        
        //- Apply collision
        virtual void collide
        (
            dsmcParcel& pP,
            dsmcParcel& pQ,
            const label cellI,
            scalar cR = -1
        );
        
        //- Apply LB energy redistribution
        virtual void redistribute
        (
            dsmcParcel& p,
            scalar& translationalEnergy,
            const scalar omegaPQ,
            const bool postReaction = false
        );
        
        //- Return the coefficients dictionary
        virtual const dictionary& coeffDict() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// #ifdef NoRepository
// #   include "LarsenBorgnakkeVariableHardSphere.C"
// #endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
