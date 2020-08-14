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
            electricalConductivityModel,
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
    electricalConductivityModel(dict, mesh),
    pe_(dict.thermo().pe()),
    localkB_(constant::physicoChemical::k.value())
{
    Info << "Loading the Spitzer-Harm electrical conductivity model" << endl;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

SpitzerHarm::~SpitzerHarm()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void SpitzerHarm::update()
{
    forAll(sigma_, cellI)
    {
        const scalar nDe = pe_[cellI]/(localkB_*T_[cellI]);
        const scalar Tpow15 = pow(T_[cellI], 1.5);
        
        if (nDe > SMALL)
        {
            sigma_[cellI] = 1.56e-4*Tpow15/log(1.23e4*Tpow15*pow(nDe, -0.5));
        }
        else
        {
            sigma_[cellI] = 0.0;
        }
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace mhd
} // End namespace Foam

// ************************************************************************* //
