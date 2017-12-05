/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2008-2009 OpenCFD Ltd.
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
    Foam::energyScalingFunction

Description

SourceFiles
    energyScalingFunction.C
    newEnergyScalingFunction.C

\*---------------------------------------------------------------------------*/

#ifndef energyScalingFunction_H
#define energyScalingFunction_H

#include "IOdictionary.H"
#include "typeInfo.H"
#include "runTimeSelectionTables.H"
#include "autoPtr.H"
#include "pairPotentialModel.H"
#include "reducedUnits.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class energyScalingFunction Declaration
\*---------------------------------------------------------------------------*/

class energyScalingFunction
{

protected:

    // Protected data

        word name_;

        dictionary energyScalingFunctionProperties_;

        const pairPotentialModel& pairPot_;

        const reducedUnits& rU_;
        

    // Private Member Functions

        //- Disallow copy construct
        energyScalingFunction(const energyScalingFunction&);

        //- Disallow default bitwise assignment
        void operator=(const energyScalingFunction&);


public:

    //- Runtime type information
    TypeName("energyScalingFunction");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            energyScalingFunction,
            dictionary,
            (
                const word& name,
                const dictionary& energyScalingFunctionProperties,
                const pairPotentialModel& pairPot,
                const reducedUnits& rU
            ),
            (name, energyScalingFunctionProperties, pairPot, rU)
        );


    // Selectors

        //- Return a reference to the selected viscosity model
        static autoPtr<energyScalingFunction> New
        (
            const word& name,
            const dictionary& energyScalingFunctionProperties,
            const pairPotentialModel& pairPot,
            const reducedUnits& rU
        );


    // Constructors

        //- Construct from components
        energyScalingFunction
        (
            const word& name,
            const dictionary& energyScalingFunctionProperties,
            const pairPotentialModel& pairPot,
            const reducedUnits& rU
        );


    // Destructor

        virtual ~energyScalingFunction()
        {}


    // Member Functions

        virtual void scaleEnergy(scalar& e, const scalar r) const = 0;

        const dictionary& energyScalingFunctionProperties() const
        {
            return energyScalingFunctionProperties_;
        }

        //- Read energyScalingFunction dictionary
        virtual bool read
        (
            const dictionary& energyScalingFunctionProperties
        ) = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
