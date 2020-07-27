/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016-2020 hyStrath
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of hyStrath, a derivative work of OpenFOAM.

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
            electricalConductivityModel,
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
    electricalConductivityModel(dict, mesh),
    n_(dict.subDict("BushCoeffs").lookupOrDefault<scalar>("n", 0.0)),
    T0_(readScalar(dict.subDict("BushCoeffs").lookup("T0"))),
    sigma0_(readScalar(dict.subDict("BushCoeffs").lookup("sigma0")))
{
    Info << "Loading the Bush electricalConductivity model" << endl;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Bush::~Bush()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Bush::update()
{
    forAll(sigma_, cellI)
    {
        sigma_[cellI] = sigma0_*pow(Tt_[cellI]/T0_, n_);
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace mhd
} // End namespace Foam

// ************************************************************************* //
