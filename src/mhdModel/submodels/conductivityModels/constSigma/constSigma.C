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

#include "constSigma.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace mhd
    {
        defineTypeNameAndDebug(constSigma, 0);

        addToRunTimeSelectionTable
        (
            conductivityModel,
            constSigma,
            mhdModel
        );

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

constSigma::constSigma
(
    const mhdModel& dict,
    const fvMesh& mesh
)
:
    conductivityModel(dict, mesh),
    sigmaValue_("sigmaValue", dimensionSet(-1,-3, 3, 0, 0, 2, 0), dict.lookup("value"))
{
    Info << "Constant conductivity value: " << sigmaValue_ << endl;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

constSigma::~constSigma()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

volScalarField constSigma::sigma() const
{
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
            sigmaValue_
           /* dimensionedScalar
            (
                "sigma",
                dimensionSet(1, -1, -3, 0, 0, 0, 0),
                sigmaValue_
            )*/
        );
    Info << "Maximum sigma: " << max(sigma);
    


    return sigma;
}


}
}

// ************************************************************************* //
