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

\*---------------------------------------------------------------------------*/

#include "energyScalingFunction.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(energyScalingFunction, 0);
defineRunTimeSelectionTable(energyScalingFunction, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

energyScalingFunction::energyScalingFunction
(
    const word& name,
    const dictionary& energyScalingFunctionProperties,
    const pairPotentialModel& pairPot,
    const reducedUnits& rU
)
:
    name_(name),
    energyScalingFunctionProperties_(energyScalingFunctionProperties),
    pairPot_(pairPot),
    rU_(rU)
    
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //


bool energyScalingFunction::read
(
    const dictionary& energyScalingFunctionProperties
)
{
    energyScalingFunctionProperties_ = energyScalingFunctionProperties;

    return true;
}

autoPtr<energyScalingFunction> energyScalingFunction::New
(
    const word& name,
    const dictionary& energyScalingFunctionProperties,
    const pairPotentialModel& pairPot,
    const reducedUnits& rU
)
{
    word energyScalingFunctionTypeName
    (
        energyScalingFunctionProperties.lookup("energyScalingFunction")
    );

    Info<< "Selecting energy scaling function "
        << energyScalingFunctionTypeName << " for "
        << name << " potential energy." << endl;

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(energyScalingFunctionTypeName);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "energyScalingFunction::New()"
        )   << "Unknown energyScalingFunction type "
            << energyScalingFunctionTypeName << nl << nl
            << "Valid  energyScalingFunctions are: " << nl
            << dictionaryConstructorTablePtr_->toc()
            << exit(FatalError);
    }

    return autoPtr<energyScalingFunction>
        (cstrIter()(name, energyScalingFunctionProperties, pairPot, rU));
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
