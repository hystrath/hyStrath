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

#include "multiSpeciesTransportModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::multiSpeciesTransportModel>
Foam::multiSpeciesTransportModel::New
(
    rho2ReactionThermo& thermo,
    const compressible::turbulenceModel& turbulence
)
{
    const word partialModelName =
        word
        (
            thermo.transportDictionary().subDict("transportModels")
                .lookup("multiSpeciesTransport")
        );

    const word modelName = partialModelName + '<'
        + thermo.partialThermoName() + '>';

    Info<< "Loading the multispecies transport model:" << tab
        << partialModelName << "\n" << endl;

    fvMeshConstructorTable::iterator cstrIter =
        fvMeshConstructorTablePtr_->find(modelName);

    if (cstrIter == fvMeshConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "DiffusionModel::New(const volVectorField&, "
            "const surfaceScalarField&)"
        )   << "Unknown diffusionModel type "
            << modelName << endl << endl
            << "Valid diffusionModels are : " << endl
            << fvMeshConstructorTablePtr_->toc()
            << exit(FatalError);
    }
    
    return autoPtr<multiSpeciesTransportModel>
        (cstrIter()(thermo, turbulence));
}


// ************************************************************************* //
