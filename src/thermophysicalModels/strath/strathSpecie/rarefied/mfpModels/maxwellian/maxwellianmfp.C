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

#include "maxwellianmfp.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace basicMfpModels
    {
        defineTypeNameAndDebug(maxwellianmfp, 0);
        addToRunTimeSelectionTable
        (
            basicMfpModel,
            maxwellianmfp, 
            dictionary
        );
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::basicMfpModels::maxwellianmfp::maxwellianmfp
(
    const word& name,
    const label& speciesIndex,
    const dictionary& dict,
    const dictionary& dictThermoPhy,
    const volScalarField& p,
    const volScalarField& Tt
)
:
    basicMfpModel(name, speciesIndex, dict, dictThermoPhy, p, Tt)
{   
    molW_ = 1.0e-3*readScalar(dictThermoPhy.subDict(name).subDict("specie").lookup("molWeight"));
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::basicMfpModels::maxwellianmfp::mfp() const
{
    const fvMesh& mesh = this->Tt_.mesh();

    tmp<volScalarField> tmfp
    (
        new volScalarField
        (
            IOobject
            (
                "mfp_" + name_,
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionSet(0, 1, 0, 0, 0)
        )
    );

    volScalarField& mfp = tmfp.ref();

    forAll(this->Tt_, celli)
    {
        mfp[celli] = sqrt(constant::mathematical::pi*molW_/(2.0*constant::physicoChemical::R.value()*this->Tt_[celli]));
    }
    

    forAll(this->Tt_.boundaryField(), patchi)
    {
        const fvPatchScalarField& pTt = this->Tt_.boundaryField()[patchi];
        fvPatchScalarField& pmfp = mfp.boundaryFieldRef()[patchi];

        forAll(pTt, facei)
        {
            pmfp[facei] = sqrt(constant::mathematical::pi*molW_/(2.0*constant::physicoChemical::R.value()*pTt[facei]));
        }
    }

    return tmfp;
}


Foam::tmp<Foam::scalarField> Foam::basicMfpModels::maxwellianmfp::mfp
(
    const label patchi,
    const scalarField& p,
    const scalarField& Tt
) const
{
    tmp<scalarField> tmfp(new scalarField(Tt.size()));
    scalarField& mfp = tmfp.ref();

    forAll(Tt, facei)
    {
        mfp[facei] = sqrt(constant::mathematical::pi*molW_/(2.0*constant::physicoChemical::R.value()*Tt[facei]));
    }

    return tmfp;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
