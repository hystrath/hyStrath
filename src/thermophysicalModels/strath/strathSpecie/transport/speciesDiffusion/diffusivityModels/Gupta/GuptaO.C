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

#include "GuptaO.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace binaryDiffusivityModels
    {
        defineTypeNameAndDebug(GuptaO, 0);
        addToRunTimeSelectionTable
        (
            binaryDiffusivityModel,
            GuptaO, 
            dictionary
        );
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::binaryDiffusivityModels::GuptaO::GuptaO
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
    binaryDiffusivityModel(name1, name2, dictThermo, dictTransport, p, pe, T),

    W1_(readScalar(dictThermo.subDict(name1).subDict("specie").lookup("molWeight"))*1.0e-3),
    W2_(readScalar(dictThermo.subDict(name2).subDict("specie").lookup("molWeight"))*1.0e-3),
    pi(Foam::constant::mathematical::pi),
    kB(Foam::constant::physicoChemical::k.value()),
    Runi(Foam::constant::physicoChemical::R.value()),
    constantFactorInCollisionTerm_(8.0e-20*sqrt(2.0*W1_*W2_/(pi*Runi*(W1_+W2_)))),
    e4OverkB2_(pow(4.8032e-10, 4.0)/sqr(kB))
{
    word year = word::null;
    
    if(dictTransport.subDict("transportModels")
        .subDict("diffusiveFluxesParameters").found("yearGuptaModel"))
    {
        year = word(dictTransport.subDict("transportModels")
            .subDict("diffusiveFluxesParameters").lookup("yearGuptaModel"));
    }
    else
    {
        FatalErrorIn("void Foam::binaryDiffusivityModels::GuptaO::GuptaO(...)")
            << "Entry 'yearGuptaModel' is missing in transportModels/diffusiveFluxesParameters."
            << exit(FatalError);
    }
    
    FixedList<scalar,4> defaultList;
    forAll(defaultList, i)
    {
        defaultList[i] = 0.0;
    }
    
    if(dictTransport.subDict("collisionData").subDict("neutralNeutralInteractions")
           .subDict("Gupta"+year+"O").subDict("Omega11").found(name1+"_"+name2))
    {
        piOmega1_ = dictTransport.subDict("collisionData").subDict("neutralNeutralInteractions")
           .subDict("Gupta"+year+"O").subDict("Omega11").lookupOrDefault<FixedList<scalar,4>>(name1+"_"+name2, defaultList);    
    }
    else if(dictTransport.subDict("collisionData").subDict("neutralNeutralInteractions")
           .subDict("Gupta"+year+"O").subDict("Omega11").found(name2+"_"+name1))
    {
        piOmega1_ = dictTransport.subDict("collisionData").subDict("neutralNeutralInteractions")
           .subDict("Gupta"+year+"O").subDict("Omega11").lookupOrDefault<FixedList<scalar,4>>(name2+"_"+name1, defaultList);    
    }
    else
    {
        FatalErrorIn("void Foam::binaryDiffusivityModels::GuptaO::GuptaO(...)")
            << "Gupta's 1990 curve fit data missing for species couple (" << name1 << ", " << name2 << ")."
            << exit(FatalError);
    }
}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::binaryDiffusivityModels::GuptaO::D() const
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
        d[celli] = kB*this->T_[celli]
            /(this->p_[celli]*collisionTerm1(this->T_[celli], this->pe_[celli]));
    }

    forAll(this->T_.boundaryField(), patchi)
    {
        const fvPatchScalarField& pT = this->T_.boundaryField()[patchi];
	      const fvPatchScalarField& pp = this->p_.boundaryField()[patchi];
	      const fvPatchScalarField& ppe = this->pe_.boundaryField()[patchi];
        fvPatchScalarField& pD = d.boundaryFieldRef()[patchi];

        forAll(pT, facei)
        {
            pD[facei] = kB*pT[facei]
                /(pp[facei]*collisionTerm1(pT[facei], ppe[facei]));
        }
    }

    return tD;
}


Foam::tmp<Foam::scalarField> Foam::binaryDiffusivityModels::GuptaO::D
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
        d[facei] = kB*T[facei]
            /(p[facei]*collisionTerm1(T[facei]));
    }

    return tD;
}


Foam::tmp<Foam::scalarField> Foam::binaryDiffusivityModels::GuptaO::D
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
        d[facei] = kB*T[facei]
            /(p[facei]*collisionTerm1(T[facei], pe[facei]));
    }

    return tD;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
