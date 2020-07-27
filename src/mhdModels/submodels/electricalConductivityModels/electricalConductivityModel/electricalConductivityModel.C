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

#include "electricalConductivityModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace mhd
    {
        defineTypeNameAndDebug(electricalConductivityModel, 0);
        defineRunTimeSelectionTable(electricalConductivityModel, mhdModel);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

electricalConductivityModel::electricalConductivityModel
(
    const mhdModel& dict,
    const fvMesh& mesh
)
:
    dict_(dict),
    mesh_(mesh),
    Tt_(dict.thermo().Tt()),
    sigma_
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
        dimensionedScalar("sigma", dimensionSet(-1, -3, 3, 0, 0, 2, 0), 0.0)
    )
{}


// * * * * * * * * * * * * * * Selector * * * * * * * * * * * * * //

autoPtr<electricalConductivityModel>
electricalConductivityModel::New
(
    const mhdModel& dict,
    const fvMesh& mesh
)
{
    const word modelType(dict.lookup("electricalConductivityModel"));

    Info<< "Selecting electrical conductivity model " << modelType << endl;

    mhdModelConstructorTable::iterator cstrIter =
        mhdModelConstructorTablePtr_->find(modelType);

    if (cstrIter == mhdModelConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "electricalConductivityModel::New(const mhdModel&, const fvMesh&)"
        )   << "Unknown electricalConductivityModel type "
            << modelType << nl << nl
            << "Valid electricalConductivityModel types are: " << nl
            << mhdModelConstructorTablePtr_->sortedToc() << exit(FatalError);
    }

    return autoPtr<electricalConductivityModel>(cstrIter()(dict, mesh));
}


// * * * * * * * * * * * * * * * * Destructor    * * * * * * * * * * * * * * //

electricalConductivityModel::~electricalConductivityModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} //End namespace mhd
} //End namespace Foam

// ************************************************************************* //
