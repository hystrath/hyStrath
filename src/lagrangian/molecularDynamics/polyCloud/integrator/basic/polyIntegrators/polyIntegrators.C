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

#include "polyIntegrators.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
//- Null Constructor (for all other md constructors)
polyIntegrators::polyIntegrators
(
    Time& t,
    const polyMesh& mesh
)
:
    time_(t),
    mesh_(mesh),
    controlDict_
    (
        IOobject
        (
            "controlDict",
            time_.system(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    )
{}

//- Constructor for mdFOAM
polyIntegrators::polyIntegrators
(
    Time& t,
    const polyMesh& mesh,
    polyMoleculeCloud& molCloud
)
:
    time_(t),
    mesh_(mesh),
    controlDict_
    (
        IOobject
        (
            "controlDict",
            time_.system(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    )
{
    
    int_ = autoPtr<polyIntegrator>
    (
        polyIntegrator::New(time_, molCloud, controlDict_)
    );
}

polyIntegrators::~polyIntegrators()
{}

autoPtr<polyIntegrator>& polyIntegrators::integrator()
{
    return int_;
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
