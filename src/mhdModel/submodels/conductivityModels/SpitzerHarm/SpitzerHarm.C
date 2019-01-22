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

#include "SpitzerHarm.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace mhd
    {
        defineTypeNameAndDebug(SpitzerHarm, 0);

        addToRunTimeSelectionTable
        (
            conductivityModel,
            SpitzerHarm,
            mhdModel
        );

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

SpitzerHarm::SpitzerHarm
(
    const mhdModel& dict,
    const fvMesh& mesh
)
:
    conductivityModel(dict, mesh)
{
    Info << "Creating SpitzerHarm conductivity model" << endl;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

SpitzerHarm::~SpitzerHarm()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

volScalarField SpitzerHarm::sigma() const
{
    volScalarField nDe = dict_.thermo().composition().nD("e-");
    volScalarField T = dict_.thermo().Tt();
    volScalarField sigma
    (
        IOobject
        (
            "sigma",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar
        (
            "sigma",
            dimensionSet(1, -1, -3, 0, 0, 0, 0),
            0.0
        )
    );
    forAll(sigma, cellI)
    {
        if(nDe[cellI] == 0.0)
        {
            sigma[cellI] = 0.0;
        }
        else
        {
            sigma[cellI] = 1.56e-4*pow(T[cellI], 1.5)/log(1.23e+4*pow(T[cellI], 1.5)*pow(nDe[cellI], -0.5));
        }
    }
    return sigma;
}


}
}

// ************************************************************************* //
