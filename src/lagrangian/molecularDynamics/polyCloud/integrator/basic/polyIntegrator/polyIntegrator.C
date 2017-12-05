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

#include "polyIntegrator.H"
#include "IFstream.H"
#include "graph.H"
#include "polyMoleculeCloud.H"


namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(polyIntegrator, 0);

defineRunTimeSelectionTable(polyIntegrator, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
polyIntegrator::polyIntegrator
(
    Time& t,
    polyMoleculeCloud& molCloud,
    const dictionary& dict
)
:
    mesh_(refCast<const fvMesh>(molCloud.mesh())),
    molCloud_(molCloud),
    time_(t)

{
}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

autoPtr<polyIntegrator> polyIntegrator::New
(
    Time& t,
    polyMoleculeCloud& molCloud,
    const dictionary& dict
)
{
    word polyIntegratorName = "velocityVerlet";
    
    if(dict.found("integrator"))
    {
        const word polyIntegratorNameTemp
        (
            dict.lookup("integrator")
        );
        
        polyIntegratorName = polyIntegratorNameTemp;
    }

    Info<< "Selecting polyIntegrator "
         << polyIntegratorName << endl;

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(polyIntegratorName);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalError
            << "polyIntegrator::New(const dictionary&) : " << endl
            << "    unknown polyIntegrator type "
            << polyIntegratorName
            << ", constructor not in hash table" << endl << endl
            << "    Valid types are :" << endl;
        Info<< dictionaryConstructorTablePtr_->toc() << abort(FatalError);
    }

    return autoPtr<polyIntegrator>
    (
        cstrIter()(t, molCloud, dict)
    );
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

polyIntegrator::~polyIntegrator()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //





} // End namespace Foam

// ************************************************************************* //
