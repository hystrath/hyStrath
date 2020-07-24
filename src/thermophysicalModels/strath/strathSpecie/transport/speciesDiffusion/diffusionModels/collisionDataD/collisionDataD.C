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
    namespace binaryDiffusionModels
    {
        defineTypeNameAndDebug(collisionDataD, 0);
        addToRunTimeSelectionTable
        (
            binaryDiffusionModel,
            collisionDataD,
            dictionary
        );
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::binaryDiffusionModels::collisionDataD::collisionDataD
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
    binaryDiffusionModel(name1, name2, dictThermo, dictTransport, p, pe, T)
{
    word collisionDataModel = word::null;

    if
    (
        dictTransport.subDict("transportModels")
            .subDict("diffusionModelParameters").found("collisionDataModel")
    )
    {
        collisionDataModel =
            word
            (
                dictTransport.subDict("transportModels")
                    .subDict("diffusionModelParameters")
                    .lookup("collisionDataModel")
            );
    }
    else
    {
        FatalErrorIn
        (
            "collisionDataD::collisionDataD(const word&, const word&, "
            "const dictionary&, const dictionary&, const volScalarField&, "
            "const volScalarField&, const volScalarField&)"
        )   << "Entry 'collisionDataModel' is missing in transportModels/"
            << "diffusionModelParameters."
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
            .subDict("tabulatedInteractions").subDict(collisionDataModel)
            .subDict("Dbar").found(name1+"_"+name2)
    )
    {
        Dbar_ = dictTransport.subDict("collisionData")
            .subDict("tabulatedInteractions")
            .subDict(collisionDataModel).subDict("Dbar")
            .lookupOrDefault<FixedList<scalar,4>>(name1+"_"+name2, defaultList);
    }
    else if
    (
        dictTransport.subDict("collisionData")
            .subDict("tabulatedInteractions").subDict(collisionDataModel)
            .subDict("Dbar").found(name2+"_"+name1)
    )
    {
        Dbar_ = dictTransport.subDict("collisionData")
            .subDict("tabulatedInteractions")
            .subDict(collisionDataModel).subDict("Dbar")
            .lookupOrDefault<FixedList<scalar,4>>(name2+"_"+name1, defaultList);
    }
    else
    {
        FatalErrorIn
        (
            "collisionDataD::collisionDataD(const word&, const word&, "
            "const dictionary&, const dictionary&, const volScalarField&, "
            "const volScalarField&, const volScalarField&)"
        )   << "Collision integral data missing for species couple (" 
            << name1 << ", " << name2 << ")."
            << exit(FatalError);
    }

    Dbar_[3] = Pstd*exp(Dbar_[3])/1.0e4;
}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::binaryDiffusionModels::collisionDataD::D() const
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
    
    const scalarField& T = this->T_.internalField();
    const scalarField& p = this->p_.internalField();
    const scalarField& pe = this->pe_.internalField();

    forAll(T, celli)
    {
        d[celli] = DijBar(T[celli], pe[celli])/p[celli];
    }

    forAll(this->T_.boundaryField(), patchi)
    {
        const fvPatchScalarField& pT = this->T_.boundaryField()[patchi];
	      const fvPatchScalarField& pp = this->p_.boundaryField()[patchi];
	      const fvPatchScalarField& ppe = this->pe_.boundaryField()[patchi];
	      
        fvPatchScalarField& pD = d.boundaryFieldRef()[patchi];

        forAll(pT, facei)
        {
            pD[facei] = DijBar(pT[facei], ppe[facei])/pp[facei];
        }
    }

    return tD;
}


Foam::tmp<Foam::scalarField>
Foam::binaryDiffusionModels::collisionDataD::D
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
