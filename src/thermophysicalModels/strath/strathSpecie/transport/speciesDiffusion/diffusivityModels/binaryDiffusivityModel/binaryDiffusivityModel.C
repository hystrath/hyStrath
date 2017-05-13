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

#include "binaryDiffusivityModel.H"
#include "volFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(binaryDiffusivityModel, 0);
    defineRunTimeSelectionTable(binaryDiffusivityModel, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::binaryDiffusivityModel::binaryDiffusivityModel
(
    const word& name1,
    const word& name2,
    const dictionary& dictThermo,
    const dictionary& dictTransport,
    const volScalarField& p,
    const volScalarField& T
)
:
    name1_(name1),    
    name2_(name2),
    dictThermo_(dictThermo),
    dictTransport_(dictTransport),    
    p_(p),
    pe_(p*0.0),  
    T_(T)
{}


Foam::binaryDiffusivityModel::binaryDiffusivityModel
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
    name1_(name1),    
    name2_(name2),
    dictThermo_(dictThermo),
    dictTransport_(dictTransport),    
    p_(p),  
    pe_(pe),  
    T_(T)
{
    if(name1_[-1] == '-')
    {
        if(name2_[-1] == '-') collisionType_ = 4;
        else if(name2_[-1] == '+') collisionType_ = 3;
        else collisionType_ = 1;
    }
    else if(name1_[-1] == '+')
    {
        if(name2_[-1] == '-') collisionType_ = 3;
        else if(name2_[-1] == '+') collisionType_ = 2;
        else collisionType_ = 1;
    }
    else
    {
        if(name2_[-1] == '-' or name2_[-1] == '+') collisionType_ = 1;
        else collisionType_ = 0;    
    }
}

// * * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::binaryDiffusivityModel> Foam::binaryDiffusivityModel::New
(
    const word& name1,
    const word& name2,
    const dictionary& dictThermo,
    const dictionary& dictTransport,
    const volScalarField& p,
    const volScalarField& T
)
{
    word binaryDiffusivityModelTypeName
         (
            dictTransport.subDict("transportModels")
                .lookup("binaryDiffusivityModel")
         );

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(binaryDiffusivityModelTypeName);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "DiffusivityModel::New(const volVectorField&, "
            "const surfaceScalarField&)"
        )   << "Unknown binaryDiffusivityModel type "
            << binaryDiffusivityModelTypeName << endl << endl
            << "Valid  binaryDiffusivityModels are : " << endl
            << dictionaryConstructorTablePtr_->toc()
            << exit(FatalError);
    }

    return autoPtr<binaryDiffusivityModel>
        (cstrIter()(name1, name2, dictThermo, dictTransport, p, p*0.0, T));
}


Foam::autoPtr<Foam::binaryDiffusivityModel> Foam::binaryDiffusivityModel::New
(
    const word& name1,
    const word& name2,
    const dictionary& dictThermo,
    const dictionary& dictTransport,
    const volScalarField& p,
    const volScalarField& pe,
    const volScalarField& T
)
{
    word binaryDiffusivityModelTypeName
         (
            dictTransport.subDict("transportModels")
                .lookup("binaryDiffusivityModel")
         );

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(binaryDiffusivityModelTypeName);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "DiffusivityModel::New(const volVectorField&, "
            "const surfaceScalarField&)"
        )   << "Unknown binaryDiffusivityModel type "
            << binaryDiffusivityModelTypeName << endl << endl
            << "Valid  binaryDiffusivityModels are : " << endl
            << dictionaryConstructorTablePtr_->toc()
            << exit(FatalError);
    }

    return autoPtr<binaryDiffusivityModel>
        (cstrIter()(name1, name2, dictThermo, dictTransport, p, pe, T));
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::binaryDiffusivityModel::D() const
{
    const fvMesh& mesh = T_.mesh();

    tmp<volScalarField> tD
    (
        new volScalarField
        (
            IOobject
            (
                "D_" + name1_ + "_" + name2_,
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionSet(0, 2, -1, 0, 0)
        )
    );

    volScalarField& d = tD();

    forAll(T_, celli)
    {
        d[celli] = Foam::VSMALL;
    }

    forAll(T_.boundaryField(), patchi)
    {
        const fvPatchScalarField& pT = T_.boundaryField()[patchi];
        fvPatchScalarField& pD = d.boundaryField()[patchi];

        forAll(pT, facei)
        {
            pD[facei] = Foam::VSMALL;
        }
    }

    return tD;
}


Foam::tmp<Foam::scalarField> Foam::binaryDiffusivityModel::D
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
        d[facei] = Foam::VSMALL;
    }

    return tD;
}


Foam::tmp<Foam::scalarField> Foam::binaryDiffusivityModel::D
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
        d[facei] = Foam::VSMALL;
    }

    return tD;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
