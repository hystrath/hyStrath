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

#include "GuptaD.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace binaryDiffusivityModels
    {
        defineTypeNameAndDebug(GuptaD, 0);
        addToRunTimeSelectionTable
        (
            binaryDiffusivityModel,
            GuptaD, 
            dictionary
        );
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::binaryDiffusivityModels::GuptaD::GuptaD
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
    word year = word::null;
    
    if(dictTransport.subDict("transportModels")
        .subDict("diffusiveFluxesParameters").found("yearGuptaModel"))
    {
        year = word(dictTransport.subDict("transportModels")
            .subDict("diffusiveFluxesParameters").lookup("yearGuptaModel"));
    }
    else
    {
        FatalErrorIn("void Foam::binaryDiffusivityModels::GuptaD::GuptaD(...)")
            << "Entry 'yearGuptaModel' is missing in transportModels/diffusiveFluxesParameters."
            << exit(FatalError);
    }
    
    if(dictTransport.subDict("collisionData").subDict("neutralNeutralInteractions")
           .subDict("Gupta"+year+"D").subDict("Dbar").found(name1+"_"+name2))
    {
        Dbar_ = dictTransport.subDict("collisionData").subDict("neutralNeutralInteractions")
           .subDict("Gupta"+year+"D").subDict("Dbar").lookup(name1+"_"+name2);    
    }
    else if(dictTransport.subDict("collisionData").subDict("neutralNeutralInteractions")
           .subDict("Gupta"+year+"D").subDict("Dbar").found(name2+"_"+name1))
    {
        Dbar_ = dictTransport.subDict("collisionData").subDict("neutralNeutralInteractions")
           .subDict("Gupta"+year+"D").subDict("Dbar").lookup(name2+"_"+name1);    
    }
    else
    {
        FatalErrorIn("void Foam::binaryDiffusivityModels::GuptaD::GuptaD(...)")
            << "Collision integral data missing for species couple (" << name1 << ", " << name2 << ")."
            << exit(FatalError);
    }
    
    Dbar_[3] = 1.01325e5*exp(Dbar_[3])/1.0e4;
}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::binaryDiffusivityModels::GuptaD::D() const
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
        d[celli] = DijBar(this->T_[celli], this->pe_[celli])/this->p_[celli];
    }

    forAll(this->T_.boundaryField(), patchi)
    {
        const fvPatchScalarField& pT = this->T_.boundaryField()[patchi];
	      const fvPatchScalarField& pp = this->p_.boundaryField()[patchi];
	      const fvPatchScalarField& ppe = this->pe_.boundaryField()[patchi];
        fvPatchScalarField& pD = d.boundaryField()[patchi];

        forAll(pT, facei)
        {
            pD[facei] = DijBar(pT[facei], ppe[facei])/pp[facei];
        }
    }

    return tD;
}


Foam::tmp<Foam::scalarField> 
Foam::binaryDiffusivityModels::GuptaD::D
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


Foam::tmp<Foam::scalarField> 
Foam::binaryDiffusivityModels::GuptaD::D
(
    const scalarField& p,
    const scalarField& pe,
    const scalarField& T,
    const label patchi
) const
{
    tmp<scalarField> tD(new scalarField(T.size()));
    scalarField& d = tD();

    forAll(T, facei)
    {
        d[facei] = DijBar(T[facei], pe[facei])/p[facei];
    }

    return tD;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
