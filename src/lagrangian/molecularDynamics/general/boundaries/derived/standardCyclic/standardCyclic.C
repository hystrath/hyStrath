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

Description

\*---------------------------------------------------------------------------*/

#include "standardCyclic.H"
#include "addToRunTimeSelectionTable.H"
#include "fvc.H"
#include "processorPolyPatch.H"
#include "cyclicPolyPatch.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(standardCyclic, 0);

addToRunTimeSelectionTable(cyclicBoundary, standardCyclic, dictionary);



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
standardCyclic::standardCyclic
(
    Time& t,
    const polyMesh& mesh,
    const dictionary& dict
)
:
    cyclicBoundary(t, mesh,  dict)
{
    if
    (
           (theta_ < (constant::mathematical::pi - SMALL))
        || (theta_ > (constant::mathematical::pi + SMALL))
    )
    {
        
        FatalErrorIn("standardCyclicBoundary::standardCyclicBoundary()")
            << "Patch: " << patchName_ << " needs a rotational cyclic boundary model. " 
             << " Angle = " << theta_
            << nl << "in: "
            << t.system()/"boundariesDict"
            << exit(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

standardCyclic::~standardCyclic()
{}



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //






} // End namespace Foam

// ************************************************************************* //
