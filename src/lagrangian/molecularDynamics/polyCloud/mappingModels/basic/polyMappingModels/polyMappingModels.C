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

Description

\*---------------------------------------------------------------------------*/

#include "polyMappingModels.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

//- Constructor
polyMappingModels::polyMappingModels
(
    const polyMesh& mesh,
    polyMoleculeCloud& molCloud
)
:
    polyMappingModelsDict_
    (
        IOobject
        (
            "mapLagrangianFieldsDict",
            mesh.time().system(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        )
    ),
	modelList_(polyMappingModelsDict_.lookup("poly")),
	models_(modelList_.size())
{
    if( models_.size() > 0 )
    {
        forAll(models_, m)
        {
            const entry& polyMappingModelsI = modelList_[m];
            const dictionary& polyMappingModelsIDict = polyMappingModelsI.dict();

            Info << nl << "Mapping polyMolecules from model #: " << m << endl;

            models_[m] = autoPtr<polyMappingModel>
            (
                polyMappingModel::New(molCloud, polyMappingModelsIDict)
            );
        }
    }
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
