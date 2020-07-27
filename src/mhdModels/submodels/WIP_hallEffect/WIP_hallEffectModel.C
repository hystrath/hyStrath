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

#include "error.H"
#include "hallEffectModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace mhd
    {
        defineTypeNameAndDebug(hallEffectModel, 0);
        defineRunTimeSelectionTable(hallEffectModel, mhdModel);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

hallEffectModel::hallEffectModel
(
    const volVectorField& B,
    const fvMesh& mesh
)
:
    B_(B),
    mesh_(mesh)
{}


// * * * * * * * * * * * * * * Selector * * * * * * * * * * * * * //

autoPtr<hallEffectModel>
hallModel::New
(
    const mhdModel& dict,
    const fvMesh& mesh
)
{
    const word modelType(dict.lookup("conductivityModel"));

    Info<< "Selecting conductivityModel " << modelType << endl;

    mhdModelConstructorTable::iterator cstrIter =
        mhdModelConstructorTablePtr_->find(modelType);

    if (cstrIter == mhdModelConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "conductivityModel::New(const mhdModel&, const fvMesh&)"
        )   << "Unknown conductivityModel type "
            << modelType << nl << nl
            << "Valid conductivityModel types are: " << nl
            << mhdModelConstructorTablePtr_->sortedToc() << exit(FatalError);
    }

    return autoPtr<conductivityModel>(cstrIter()(dict, mesh));
}


// * * * * * * * * * * * * * * * * Destructor    * * * * * * * * * * * * * * //

conductivityModel::~conductivityModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

volScalarField conductivityModel::sigma() const
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
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE,
                false
            ),
            mesh_,
            dimensionedScalar("zero", dimless/dimLength, 0.0)
        )
    );
}


void conductivityModel::correct()
{

}


} //End namespace mhd
} //End namespace Foam

// ************************************************************************* //
