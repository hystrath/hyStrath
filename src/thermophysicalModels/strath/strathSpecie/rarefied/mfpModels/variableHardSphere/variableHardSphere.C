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

#include "variableHardSphere.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace basicMfpModels
    {
        defineTypeNameAndDebug(variableHardSphere, 0);
        addToRunTimeSelectionTable
        (
            basicMfpModel,
            variableHardSphere,
            dictionary
        );
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::basicMfpModels::variableHardSphere::variableHardSphere
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
    omega_ = readScalar(dictThermoPhy.subDict(name).subDict("specie").lookup("omega"));
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::basicMfpModels::variableHardSphere::mfp() const
{
    const fvMesh& mesh = this->T_.mesh();

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

    forAll(this->T_, celli)
    {
        mfp[celli] = 2.0/15.0*(5.0-2.0*omega_)*(7.0-2.0*omega_)*sqrt(molW_/(constant::mathematical::twoPi*constant::physicoChemical::R.value()*this->T_[celli]));
    }


    forAll(this->T_.boundaryField(), patchi)
    {
        const fvPatchScalarField& pTt = this->T_.boundaryField()[patchi];
        fvPatchScalarField& pmfp = mfp.boundaryFieldRef()[patchi];

        forAll(pTt, facei)
        {
            pmfp[facei] = 2.0/15.0*(5.0-2.0*omega_)*(7.0-2.0*omega_)*sqrt(molW_/(constant::mathematical::twoPi*constant::physicoChemical::R.value()*pTt[facei]));
        }
    }

    return tmfp;
}


Foam::tmp<Foam::scalarField> Foam::basicMfpModels::variableHardSphere::mfp
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
        mfp[facei] = 2.0/15.0*(5.0-2.0*omega_)*(7.0-2.0*omega_)*sqrt(molW_/(constant::mathematical::twoPi*constant::physicoChemical::R.value()*Tt[facei]));
    }

    return tmfp;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
