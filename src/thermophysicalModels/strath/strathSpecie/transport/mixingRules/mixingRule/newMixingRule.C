/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
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

#include "mixingRule.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::mixingRule>
Foam::mixingRule::New
(
    rho2ReactionThermo& thermo,
    const compressibleTurbulenceModel& turbulence
)
{
    const word partialMixingRuleName = word(thermo.transportDictionary()
            .subDict("transportModels").lookup("mixingRule"));
        
    word mixingRuleName = partialMixingRuleName + word("MR") 
        +'<' + thermo.partialThermoName() + '>'; 
        
    Info<< "\nLoading the transport mixing rule:" << tab << partialMixingRuleName 
        << "\n" << endl;

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
