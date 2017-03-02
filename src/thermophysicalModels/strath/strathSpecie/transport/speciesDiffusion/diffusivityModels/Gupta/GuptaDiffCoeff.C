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

#include "GuptaDiffCoeff.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace binaryDiffusivityModels
    {
        defineTypeNameAndDebug(GuptaDiffCoeff, 0);
        addToRunTimeSelectionTable
        (
            binaryDiffusivityModel,
            GuptaDiffCoeff, 
            dictionary
        );
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::binaryDiffusivityModels::GuptaDiffCoeff::GuptaDiffCoeff
(
    const word& name1,
    const word& name2,
    const dictionary& dictThermo,
    const dictionary& dictTransport,
    const volScalarField& p,
    const volScalarField& T
)
:
    binaryDiffusivityModel(name1, name2, dictThermo, dictTransport, p, T),

    W1_(readScalar(dictThermo.subDict(name1).subDict("specie").lookup("molWeight"))*1.0e-3),
    W2_(readScalar(dictThermo.subDict(name2).subDict("specie").lookup("molWeight"))*1.0e-3),
    pi(Foam::constant::mathematical::pi),
    kB(Foam::constant::physicoChemical::k.value()),
    Runi(Foam::constant::physicoChemical::R.value())
{
    if (dictTransport.subDict("collisionIntegrals").subDict("involvingNeutral")
           .subDict("GuptaDiffCoeff_Omega11").found(name1+"_"+name2))
    {
        piOmega_ = dictTransport.subDict("collisionIntegrals").subDict("involvingNeutral")
           .subDict("GuptaDiffCoeff_Omega11").lookup(name1+"_"+name2);    
    }
    else if (dictTransport.subDict("collisionIntegrals").subDict("involvingNeutral")
            .subDict("GuptaDiffCoeff_Omega11").found(name2+"_"+name1))
    {
        piOmega_ = dictTransport.subDict("collisionIntegrals").subDict("involvingNeutral")
            .subDict("GuptaDiffCoeff_Omega11").lookup(name2+"_"+name1);    
    }
    else
    {
        FatalErrorIn("void Foam::binaryDiffusivityModels::GuptaDiffCoeff::GuptaDiffCoeff(...)")
            << "Collision integral data missing for species couple (" << name1 << ", " << name2 << ")."
            << exit(FatalError);
    }
}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::binaryDiffusivityModels::GuptaDiffCoeff::D() const
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

    volScalarField& d = tD();

    forAll(this->T_, celli)
    {
        d[celli] = DijBar(this->T_[celli])/this->p_[celli];
    }

    forAll(this->T_.boundaryField(), patchi)
    {
        const fvPatchScalarField& pT = this->T_.boundaryField()[patchi];
	      const fvPatchScalarField& pp = this->p_.boundaryField()[patchi];
        fvPatchScalarField& pD = d.boundaryField()[patchi];

        forAll(pT, facei)
        {
            pD[facei] = DijBar(pT[facei])/pp[facei];
        }
    }

    return tD;
}


Foam::tmp<Foam::scalarField> Foam::binaryDiffusivityModels::GuptaDiffCoeff::D
(
    const scalarField& p,
    const scalarField& T,
    const label patchi
) const
{
    tmp<scalarField> tD(new scalarField(T.size()));
    scalarField& d = tD();

    forAll(T, facei)
    {
        d[facei] = DijBar(T[facei])/p[facei];
    }

    return tD;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
