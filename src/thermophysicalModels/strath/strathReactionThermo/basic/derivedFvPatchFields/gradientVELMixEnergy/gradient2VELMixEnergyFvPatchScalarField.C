/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
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

#include "gradient2VELMixEnergyFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "multi2Thermo.H" // NEW VINCENT

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::gradient2VELMixEnergyFvPatchScalarField::
gradient2VELMixEnergyFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(p, iF)
{} // Only this constructor is used


Foam::gradient2VELMixEnergyFvPatchScalarField::
gradient2VELMixEnergyFvPatchScalarField
(
    const gradient2VELMixEnergyFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedGradientFvPatchScalarField(ptf, p, iF, mapper)
{}


Foam::gradient2VELMixEnergyFvPatchScalarField::
gradient2VELMixEnergyFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedGradientFvPatchScalarField(p, iF, dict)
{}


Foam::gradient2VELMixEnergyFvPatchScalarField::
gradient2VELMixEnergyFvPatchScalarField
(
    const gradient2VELMixEnergyFvPatchScalarField& tppsf
)
:
    fixedGradientFvPatchScalarField(tppsf)
{}


Foam::gradient2VELMixEnergyFvPatchScalarField::
gradient2VELMixEnergyFvPatchScalarField
(
    const gradient2VELMixEnergyFvPatchScalarField& tppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(tppsf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::gradient2VELMixEnergyFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }
    
    //Info << "gradient2VELMixEnergy is used for patch called " << patch().name() << endl;

    const multi2Thermo& thermo = multi2Thermo::lookup2Thermo(*this);
    const label patchi = patch().index();

    const scalarField& pw = thermo.p().boundaryField()[patchi];
    
    fvPatchScalarField& Tvw =
        const_cast<fvPatchScalarField&>(thermo.Tv().boundaryField()[patchi]);
    Tvw.evaluate();
    
    gradient() = thermo.Cv_v(pw, Tvw, patchi)*Tvw.snGrad()
      + patch().deltaCoeffs()*
        (
            thermo.hevel(pw, Tvw, patchi)
          - thermo.hevel(pw, Tvw, patch().faceCells())
        );
        
    fixedGradientFvPatchScalarField::updateCoeffs();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        gradient2VELMixEnergyFvPatchScalarField
    );
}

// ************************************************************************* //
