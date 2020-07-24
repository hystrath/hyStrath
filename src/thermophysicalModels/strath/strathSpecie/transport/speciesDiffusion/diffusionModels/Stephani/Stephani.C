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

#include "Stephani.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace binaryDiffusionModels
    {
        defineTypeNameAndDebug(Stephani, 0);
        addToRunTimeSelectionTable
        (
            binaryDiffusionModel,
            Stephani,
            dictionary
        );
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::binaryDiffusionModels::Stephani::Stephani
(
    const word& name1,
    const word& name2,
    const dictionary& dictThermo,
    const dictionary& dictTransport,
    const volScalarField& p,
    const volScalarField& pe,
    const volScalarField& T
)
:
    binaryDiffusionModel(name1, name2, dictThermo, dictTransport, p, pe, T),

    pi(Foam::constant::mathematical::pi),
    kB(Foam::constant::physicoChemical::k.value()),
    Runi(Foam::constant::physicoChemical::R.value()),

    MStar_
    (
        1e-3
      * readScalar
        (
            dictThermo.subDict(name1).subDict("specie").lookup("molWeight")
        )
      * readScalar
        (
            dictThermo.subDict(name2).subDict("specie").lookup("molWeight")
        )
      / (
            readScalar
            (
                dictThermo.subDict(name1).subDict("specie").lookup("molWeight")
            ) 
          + readScalar
            (
                dictThermo.subDict(name2).subDict("specie").lookup("molWeight")
            )
        )
    ),
    mStar_(MStar_*kB/Runi),
    omegaref_
    (
        0.5*
        (
            readScalar
            (
                dictThermo.subDict(name1).subDict("specie").lookup("omega")
            )
          + readScalar
            (
                dictThermo.subDict(name2).subDict("specie").lookup("omega")
            )
        )
    ),
    dref_
    (
        0.5*
        (
            readScalar
            (
                dictThermo.subDict(name1).subDict("specie").lookup("diameter")
            )
          + readScalar
            (
                dictThermo.subDict(name2).subDict("specie").lookup("diameter")
            )
        )
    ),
    Tref_(273.0),
    constantOmegaNeutral1_(constantExpressionInOmegaNeutral1())
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::binaryDiffusionModels::Stephani::D() const
{
    const fvMesh& mesh = this->T_.mesh();

    tmp<volScalarField> tD
    (
        new volScalarField
        (
            IOobject
            (
                "rhoD_" + name1_ + "_" + name2_,
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimArea/dimTime
        )
    );

    volScalarField& d = tD.ref();

    forAll(this->T_, celli)
    {
        d[celli] = 3.0*sqr(kB*this->T_[celli])
            /(16.0*mStar_*this->p_[celli]*omegaNeutral1(this->T_[celli]));
    }

    forAll(this->T_.boundaryField(), patchi)
    {
        const fvPatchScalarField& pT = this->T_.boundaryField()[patchi];
	      const fvPatchScalarField& pp = this->p_.boundaryField()[patchi];
        fvPatchScalarField& pD = d.boundaryFieldRef()[patchi];

        forAll(pT, facei)
        {
            pD[facei] = 3.0*sqr(kB*pT[facei])
                /(16.0*mStar_*pp[facei]*omegaNeutral1(pT[facei]));
        }
    }

    return tD;
}


Foam::tmp<Foam::scalarField> Foam::binaryDiffusionModels::Stephani::D
(
    const scalarField& p,
    const scalarField& pe,
    const scalarField& T,
    const label patchi
) const
{
    tmp<scalarField> tD(new scalarField(T.size()));
    scalarField& d = tD.ref();

    forAll(T, facei)
    {
        d[facei] = 3.0*sqr(kB*T[facei])
            /(16.0*mStar_*p[facei]*omegaNeutral1(T[facei]));
    }

    return tD;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
