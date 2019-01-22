/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "Bush.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace mhd
    {
        defineTypeNameAndDebug(Bush, 0);

        addToRunTimeSelectionTable
        (
            conductivityModel,
            Bush,
            mhdModel
        );

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Bush::Bush
(
    const mhdModel& dict,
    const fvMesh& mesh
)
:
    conductivityModel(dict, mesh),
    n_(dict.lookupOrDefault<label>("n", 0)),
    T0_("T0", dimensionSet(0, 0, 0, 1, 0, 0, 0), dict.lookup("T0")),
    sigma0_("sigma0", dimensionSet(-1, -3, 3, 0, 0, 2, 0), dict.lookup("sigma0"))
{
    Info << "Creating Bush conductivity model" << endl;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Bush::~Bush()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

volScalarField Bush::sigma() const
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "sigma",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE,
                false
            ),
            sigma0_*pow(dict_.thermo().Tt()/T0_, n_)
        )
    );
}


}
}

// ************************************************************************* //
