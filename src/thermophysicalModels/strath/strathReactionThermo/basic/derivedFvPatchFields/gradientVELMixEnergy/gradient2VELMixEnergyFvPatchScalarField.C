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
    fixedGradientFvPatchScalarField(p, iF),
    thermo_(rho2ReactionThermo::lookup2ReactionThermo(*this))
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
    fixedGradientFvPatchScalarField(ptf, p, iF, mapper),
    thermo_(rho2ReactionThermo::lookup2ReactionThermo(*this))
{}


Foam::gradient2VELMixEnergyFvPatchScalarField::
gradient2VELMixEnergyFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedGradientFvPatchScalarField(p, iF, dict),
    thermo_(rho2ReactionThermo::lookup2ReactionThermo(*this))
{}


Foam::gradient2VELMixEnergyFvPatchScalarField::
gradient2VELMixEnergyFvPatchScalarField
(
    const gradient2VELMixEnergyFvPatchScalarField& tppsf
)
:
    fixedGradientFvPatchScalarField(tppsf),
    thermo_(rho2ReactionThermo::lookup2ReactionThermo(*this))
{}


Foam::gradient2VELMixEnergyFvPatchScalarField::
gradient2VELMixEnergyFvPatchScalarField
(
    const gradient2VELMixEnergyFvPatchScalarField& tppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(tppsf, iF),
    thermo_(rho2ReactionThermo::lookup2ReactionThermo(*this))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::gradient2VELMixEnergyFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }
    
    //Info << "gradient2VELMixEnergy is used for patch called " << patch().name() << endl;

    const multi2Thermo& multiThermo = multi2Thermo::lookup2Thermo(*this);
    const label patchi = patch().index();

    const scalarField& pw = multiThermo.p().boundaryField()[patchi];
    
    fvPatchScalarField& Tvw =
        const_cast<fvPatchScalarField&>(multiThermo.Tv().boundaryField()[patchi]);
    Tvw.evaluate();
    
    // NEW VINCENT 15/02/2017 *************************************************
    tmp<Field<scalar>> thevel(new Field<scalar>(pw.size()));
    Field<scalar>& hevel = thevel.ref();
    hevel = 0.0;
    
    tmp<Field<scalar>> thevelfC(new Field<scalar>(pw.size()));
    Field<scalar>& hevelfC = thevelfC.ref();
    hevelfC = 0.0;
    
    tmp<Field<scalar>> tCvvelTvw(new Field<scalar>(pw.size()));
    Field<scalar>& cvvelTvw = tCvvelTvw.ref();
    cvvelTvw = 0.0;
    
    for(label speciei=0 ; speciei<thermo_.composition().Y().size() ; speciei++)
    {
        fvPatchScalarField& spYw =
            const_cast<fvPatchScalarField&>(thermo_.composition().Y(speciei).boundaryField()[patchi]);
        spYw.evaluate();
        
        hevel += spYw*thermo_.composition().hevel(speciei, pw, Tvw, patchi);
        hevelfC += spYw*thermo_.composition().hevel(speciei, pw, Tvw, patch().faceCells());
        cvvelTvw += spYw*thermo_.composition().Cv_vel(speciei, pw, Tvw, patchi)*Tvw.snGrad();
    }
    
    gradient() = tCvvelTvw
      + patch().deltaCoeffs()*
        (
            thevel - thevelfC
        );
    // END NEW VINCENT 15/02/2017 *********************************************
        
    /*gradient() = thermo.Cv_v(pw, Tvw, patchi)*Tvw.snGrad()
      + patch().deltaCoeffs()*
        (
            thermo.hevel(pw, Tvw, patchi)
          - thermo.hevel(pw, Tvw, patch().faceCells())
        );*/ // DELETED VINCENT 15/02/2017 OLD FORMULATION
        
        
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
