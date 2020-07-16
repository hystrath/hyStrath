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

#include "mixingRule.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::mixingRule>
Foam::mixingRule::New
(
    rho2ReactionThermo& thermo,
    const compressibleTurbulenceModel& turbulence
)
{
    word partialMixingRuleName = word::null;
    
    if (thermo.composition().species().size() == 1)
    {
        partialMixingRuleName = "molar";
    }
    else
    {
        partialMixingRuleName = thermo.transportDictionary()
            .subDict("transportModels")
            .lookupOrDefault<word>("mixingRule", "Wilke");
    }
    
    word mixingRuleName = partialMixingRuleName + word("MR")
        +'<' + thermo.partialThermoName() + '>';

    Info<< "\nLoading the transport mixing rule:" << tab 
        << partialMixingRuleName << "\n" << endl;

    fvMeshConstructorTable::iterator cstrIter =
        fvMeshConstructorTablePtr_->find(mixingRuleName);

    if (cstrIter == fvMeshConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "mixingRuleModel::New(const volVectorField&, "
            "const surfaceScalarField&)"
        )   << "Unknown mixingRuleModel type "
            << mixingRuleName << endl << endl
            << "Valid  mixingRuleModels are : " << endl
            << fvMeshConstructorTablePtr_->toc()
            << exit(FatalError);
    }

    return autoPtr<mixingRule>
        (cstrIter()(thermo, turbulence));
}


// ************************************************************************* //
