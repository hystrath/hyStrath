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

#include "polyAllConfigurations.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
//- Null Constructor 
polyAllConfigurations::polyAllConfigurations
(
    const polyMesh& mesh
)
:
    mdInitialiseDict_
    (
        IOobject
        (
            "mdInitialiseDict",
            mesh.time().system(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        )
    ),

    configurationList_(),
    ids_(),
    configurations_()
{}

//- Constructor for mdInitialise
polyAllConfigurations::polyAllConfigurations
(
    const polyMesh& mesh,
    polyMoleculeCloud& molCloud
)
:
    mdInitialiseDict_
    (
        IOobject
        (
            "mdInitialiseDict",
            mesh.time().system(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        )
    ),
    configurationList_(mdInitialiseDict_.lookup("polyConfigurations")),
    ids_(configurationList_.size()),
    configurations_(configurationList_.size())
{

    Info << nl << "Creating polyConfigurations: " << nl << endl;

    if(configurations_.size() > 0 )
    {
        forAll(configurations_, c)
        {
            const entry& configurationI = configurationList_[c];
            const dictionary& configurationIDict = configurationI.dict();
    
            configurations_[c] = autoPtr<polyConfiguration>
            (
                polyConfiguration::New(molCloud, configurationIDict)
            );

            ids_[c] = c;
        }
    }
}

//- initial configuration

void polyAllConfigurations::setInitialConfig()
{
    forAll(configurations_, c)
    {
        configurations_[c]->setInitialConfiguration();
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
