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

#include "ChapmanCowling.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace mhd
    {
        defineTypeNameAndDebug(ChapmanCowling, 0);

        addToRunTimeSelectionTable
        (
            electricalConductivityModel,
            ChapmanCowling,
            mhdModel
        );
        

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

ChapmanCowling::ChapmanCowling
(
    const mhdModel& dict,
    const fvMesh& mesh
)
:
    electricalConductivityModel(dict, mesh),
    pe_(dict.thermo().pe()),
    localkB_(constant::physicoChemical::k.value())
{
    Info << "Loading the Chapman-Cowling electricalConductivity model" << endl;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

ChapmanCowling::~ChapmanCowling()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void ChapmanCowling::update()
{
    forAll(sigma_, cellI)
    {
        const scalar nDe = pe_[cellI]/(localkB_*T_[cellI]);
        sigma_[cellI] = 4.0227904e-18*nDe/sqrt(T_[cellI]);
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace mhd
} // End namespace Foam

// ************************************************************************* //
