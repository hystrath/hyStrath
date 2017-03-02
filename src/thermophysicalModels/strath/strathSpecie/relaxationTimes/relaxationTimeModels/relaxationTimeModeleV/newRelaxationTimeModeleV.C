/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

#include "relaxationTimeModeleV.H"
#include "IFstream.H" // NEW VINCENT 10/08/2016

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::relaxationTimeModeleV>
Foam::relaxationTimeModeleV::New
(
    rho2ReactionThermo& thermo,
    const compressible::turbulenceModel2& turbulence
)
{
    word eVModelName;

    // Enclose the creation of the dictionary to ensure it is deleted
    // before the relaxationTimeModeleV is created otherwise the
    // dictionary is entered in the database twice
    {
        const dictionary& thermo2TModel =
        (
            IFstream
            (
                fileName(thermo.lookup("twoTemperatureDictFile")).expand()
            )()
        );

        word partialeVModelName = word(thermo2TModel.subDict("thermalRelaxationModels").subDict("eV").lookup("relaxationType"));
        
        eVModelName = partialeVModelName +'<' + thermo.partialThermoName() + '>';
        
        Info<< "Loading the e-v relaxation time model:" << tab << partialeVModelName << "\n" << endl;
    }

    fvMeshConstructorTable::iterator cstrIter =
        fvMeshConstructorTablePtr_->find(eVModelName);

    if (cstrIter == fvMeshConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "eVModel::New(const volVectorField&, "
            "const surfaceScalarField&)"
        )   << "Unknown eVModel type "
            << eVModelName << endl << endl
            << "Valid  eVModels are : " << endl
            << fvMeshConstructorTablePtr_->toc()
            << exit(FatalError);
    }

    return autoPtr<relaxationTimeModeleV>
        (cstrIter()(thermo, turbulence));
}


// ************************************************************************* //
