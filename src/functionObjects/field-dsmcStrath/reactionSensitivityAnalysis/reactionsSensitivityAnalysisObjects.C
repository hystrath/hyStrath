/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

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

#include "reactionsSensitivityAnalysis.H"
#include "psiChemistryModel.H"
#include "rhoChemistryModel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
// Psi-based chemistry
typedef functionObjects::reactionsSensitivityAnalysis<psiChemistryModel>
    psiReactionsSensitivityAnalysisFunctionObject;

defineTemplateTypeNameAndDebugWithName
(
    psiReactionsSensitivityAnalysisFunctionObject,
    "psiReactionsSensitivityAnalysis",
    0
);

// Rho-based chemistry
typedef functionObjects::reactionsSensitivityAnalysis<rhoChemistryModel>
    rhoReactionsSensitivityAnalysisFunctionObject;

defineTemplateTypeNameAndDebugWithName
(
    rhoReactionsSensitivityAnalysisFunctionObject,
    "rhoReactionsSensitivityAnalysis",
    0
);

namespace functionObjects
{
    addToRunTimeSelectionTable
    (
        functionObject,
        psiReactionsSensitivityAnalysisFunctionObject,
        dictionary
    );

    addToRunTimeSelectionTable
    (
        functionObject,
        rhoReactionsSensitivityAnalysisFunctionObject,
        dictionary
    );
}
}


// ************************************************************************* //
