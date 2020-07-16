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

#include "relaxationTimeModelHE.H"
#include "IFstream.H" // NEW VINCENT 10/08/2016

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::relaxationTimeModelHE>
Foam::relaxationTimeModelHE::New
(
    rho2ReactionThermo& thermo,
    const compressibleTurbulenceModel& turbulence
)
{
    word HEModelName;

    // Enclose the creation of the dictionary to ensure it is deleted
    // before the relaxationTimeModelHE is created otherwise the
    // dictionary is entered in the database twice
    {
        const dictionary& thermo2TModel =
        (
            IFstream
            (
                fileName(thermo.lookup("twoTemperatureDictFile")).expand()
            )()
        );

        const word partialHEModelName = word(thermo2TModel.subDict("thermalRelaxationModels").subDict("he").lookup("relaxationType"));

        HEModelName = partialHEModelName +'<' + thermo.partialThermoName() + '>';

        Info<< "Loading the h-e relaxation time model" << tab << partialHEModelName << "\n" << endl;
    }

    fvMeshConstructorTable::iterator cstrIter =
        fvMeshConstructorTablePtr_->find(HEModelName);

    if (cstrIter == fvMeshConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "HEModel::New(const volVectorField&, "
            "const surfaceScalarField&)"
        )   << "Unknown HEModel type "
            << HEModelName << endl << endl
            << "Valid  HEModels are : " << endl
            << fvMeshConstructorTablePtr_->toc()
            << exit(FatalError);
    }


    return autoPtr<relaxationTimeModelHE>
        (cstrIter()(thermo, turbulence));
}


// ************************************************************************* //
