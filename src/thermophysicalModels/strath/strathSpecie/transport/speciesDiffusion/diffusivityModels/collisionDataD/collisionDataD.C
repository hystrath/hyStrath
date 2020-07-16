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

#include "collisionDataD.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace binaryDiffusivityModels
    {
        defineTypeNameAndDebug(collisionDataD, 0);
        addToRunTimeSelectionTable
        (
            binaryDiffusivityModel,
            collisionDataD,
            dictionary
        );
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::binaryDiffusivityModels::collisionDataD::collisionDataD
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
    binaryDiffusivityModel(name1, name2, dictThermo, dictTransport, p, pe, T)
{
    word collisionDataModel = word::null;

    if (dictTransport.subDict("transportModels")
        .subDict("diffusiveFluxesParameters").found("collisionDataModel"))
    {
        collisionDataModel = word(dictTransport.subDict("transportModels")
            .subDict("diffusiveFluxesParameters").lookup("collisionDataModel"));
    }
    else
    {
        FatalErrorIn("void Foam::binaryDiffusivityModels::collisionDataD::collisionDataD(...)")
            << "Entry 'collisionDataModel' is missing in transportModels/"
            << "diffusiveFluxesParameters."
            << exit(FatalError);
    }

    FixedList<scalar,4> defaultList;
    forAll(defaultList, i)
    {
        defaultList[i] = 0.0;
    }

    if
    (
        dictTransport.subDict("collisionData")
            .subDict("neutralNeutralInteractions").subDict(collisionDataModel)
            .subDict("Dbar").found(name1+"_"+name2)
    )
    {
        Dbar_ = dictTransport.subDict("collisionData")
            .subDict("neutralNeutralInteractions")
            .subDict(collisionDataModel).subDict("Dbar")
            .lookupOrDefault<FixedList<scalar,4>>(name1+"_"+name2, defaultList);
    }
    else if
    (
        dictTransport.subDict("collisionData")
            .subDict("neutralNeutralInteractions").subDict(collisionDataModel)
            .subDict("Dbar").found(name2+"_"+name1)
    )
    {
        Dbar_ = dictTransport.subDict("collisionData")
            .subDict("neutralNeutralInteractions")
            .subDict(collisionDataModel).subDict("Dbar")
            .lookupOrDefault<FixedList<scalar,4>>(name2+"_"+name1, defaultList);
    }
    else
    {
        FatalErrorIn("void Foam::binaryDiffusivityModels::collisionDataD::collisionDataD(...)")
            << "Collision integral data missing for species couple (" 
            << name1 << ", " << name2 << ")."
            << exit(FatalError);
    }

    Dbar_[3] = 1.01325e5*exp(Dbar_[3])/1.0e4;
}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::binaryDiffusivityModels::collisionDataD::D() const
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
        d[celli] = DijBar(this->T_[celli], 0/*this->pe_[celli]*/)/this->p_[celli];
    }

    forAll(this->T_.boundaryField(), patchi)
    {
        const fvPatchScalarField& pT = this->T_.boundaryField()[patchi];
	      const fvPatchScalarField& pp = this->p_.boundaryField()[patchi];
	      //const fvPatchScalarField& ppe = this->pe_.boundaryField()[patchi];
        fvPatchScalarField& pD = d.boundaryFieldRef()[patchi];

        forAll(pT, facei)
        {
            pD[facei] = DijBar(pT[facei], 0/*ppe[facei]*/)/pp[facei];
        }
    }

    return tD;
}


Foam::tmp<Foam::scalarField>
Foam::binaryDiffusivityModels::collisionDataD::D
(
    const scalarField& p,
    const scalarField& T,
    const label patchi
) const
{
    tmp<scalarField> tD(new scalarField(T.size()));
    scalarField& d = tD.ref();

    forAll(T, facei)
    {
        d[facei] = DijBar(T[facei])/p[facei];
    }

    return tD;
}


Foam::tmp<Foam::scalarField>
Foam::binaryDiffusivityModels::collisionDataD::D
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
        d[facei] = DijBar(T[facei], pe[facei])/p[facei];
    }

    return tD;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
