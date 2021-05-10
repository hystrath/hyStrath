/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016-2021 hyStrath
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

#include "collisionDataO.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace binaryDiffusionModels
    {
        defineTypeNameAndDebug(collisionDataO, 0);
        addToRunTimeSelectionTable
        (
            binaryDiffusionModel,
            collisionDataO,
            dictionary
        );
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::binaryDiffusionModels::collisionDataO::collisionDataO
(
    const word& name1,
    const word& name2,
    const dictionary& dictThermo,
    const dictionary& dictTransport,
    const volScalarField& p,
    const volScalarField& pe,
    const volScalarField& T,
    const volScalarField& Te
)
:
    binaryDiffusionModel(name1, name2, dictThermo, dictTransport, p, pe, T, Te),
    isTabulatedData_(true),
    localkB_(Foam::constant::physicoChemical::k.value()),
    localpi_(Foam::constant::mathematical::pi),
    e2TokB_(sqr(constant::electromagnetic::e.value()*2.998e9)/(localkB_*1.0e7))
{
    scalar W1 =
        readScalar
        (
            dictThermo.subDict(name1).subDict("specie").lookup("molWeight")
        )*1.0e-3;
        
    scalar W2 =
        readScalar
        (
            dictThermo.subDict(name2).subDict("specie").lookup("molWeight")
        )*1.0e-3;
        
    constantFactor_ = 8.0e-20/3.0
        *sqrt
        (
            2.0*W1*W2
          / (
                Foam::constant::mathematical::pi
              * Foam::constant::physicoChemical::R.value()
              * (W1+W2)
            )
        );
    
    word collisionDataModel = word::null;
    word potentialType = word::null;

    if
    (
        dictTransport.subDict("transportModels")
            .subDict("diffusionModelParameters").found("collisionDataModel")
    )
    {
        collisionDataModel = word(dictTransport.subDict("transportModels")
            .subDict("diffusionModelParameters")
            .lookup("collisionDataModel"));
    }
    else
    {
        FatalErrorIn
        (
            "collisionDataO::collisionDataO(const word&, const word&, "
            "const dictionary&, const dictionary&, const volScalarField&, "
            "const volScalarField&, const volScalarField&)"
        )   << "Entry 'collisionDataModel' is missing in transportModels/"
            << "diffusionModelParameters."
            << exit(FatalError);
    }
    
    if (collisionDataModel.find("1989") != std::string::npos)
    {
        limitingElectronPressureModel_ = 0;
    }
    else
    {
        limitingElectronPressureModel_ = 1;
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
            .subDict("Omega11").found(name1+"_"+name2)
    )
    {
        piOmega1_ = dictTransport.subDict("collisionData")
            .subDict("tabulatedInteractions")
            .subDict(collisionDataModel).subDict("Omega11")
            .lookupOrDefault<FixedList<scalar,4>>
             (
                name1+"_"+name2,
                defaultList
             );
    }
    else if
    (
        dictTransport.subDict("collisionData")
            .subDict("tabulatedInteractions")
            .subDict(collisionDataModel)
            .subDict("Omega11").found(name2+"_"+name1)
    )
    {
        piOmega1_ = dictTransport.subDict("collisionData")
            .subDict("tabulatedInteractions")
            .subDict(collisionDataModel).subDict("Omega11")
            .lookupOrDefault<FixedList<scalar,4>>
             (
                  name2+"_"+name1,
                  defaultList
             );
    }
    else
    {
        if
        (
            (collisionType_ >= 2)
        and dictTransport.subDict("collisionData")
              .subDict("tabulatedInteractions")
              .subDict(collisionDataModel)
              .subDict("Omega11").found("shieldedCoulombPotential")
        )
        {
            isTabulatedData_ = false;
            
            if (collisionType_ == 3)
            {
                potentialType = "attractive";
            }
            else
            {
                potentialType = "repulsive";
            }
            
            cCDCoeffs1_[0] = dictTransport.subDict("collisionData")
                .subDict("tabulatedInteractions")
                .subDict(collisionDataModel).subDict("Omega11")
                .subDict("shieldedCoulombPotential").subDict(potentialType)
                .lookupOrDefault<scalar>("c1", 0.0);
            cCDCoeffs1_[1] = dictTransport.subDict("collisionData")
                .subDict("tabulatedInteractions")
                .subDict(collisionDataModel).subDict("Omega11")
                .subDict("shieldedCoulombPotential").subDict(potentialType)
                .lookupOrDefault<scalar>("C1", 0.0);
            cCDCoeffs1_[2] = dictTransport.subDict("collisionData")
                .subDict("tabulatedInteractions")
                .subDict(collisionDataModel).subDict("Omega11")
                .subDict("shieldedCoulombPotential").subDict(potentialType)
                .lookupOrDefault<scalar>("D1", 0.0);
        }
        else
        {
            FatalErrorIn
            (
                "collisionDataO::collisionDataO(const word&, const word&, "
                "const dictionary&, const dictionary&, const volScalarField&, "
                "const volScalarField&, const volScalarField&)"
            )   << collisionDataModel << " curve fit data missing for species "
                << "couple (" << name1 << ", " << name2 << ")."
                << exit(FatalError);
        }
    }
    
    if
    (
        dictTransport.subDict("collisionData")
            .subDict("tabulatedInteractions").subDict(collisionDataModel)
            .subDict("Omega22").found(name1+"_"+name2)
    )
    {
        piOmega2_ = dictTransport.subDict("collisionData")
            .subDict("tabulatedInteractions")
            .subDict(collisionDataModel).subDict("Omega22")
            .lookupOrDefault<FixedList<scalar,4>>
             (
                name1+"_"+name2,
                defaultList
             );
    }
    else if
    (
        dictTransport.subDict("collisionData")
            .subDict("tabulatedInteractions")
            .subDict(collisionDataModel)
            .subDict("Omega22").found(name2+"_"+name1)
    )
    {
        piOmega2_ = dictTransport.subDict("collisionData")
            .subDict("tabulatedInteractions")
            .subDict(collisionDataModel).subDict("Omega22")
            .lookupOrDefault<FixedList<scalar,4>>
             (
                  name2+"_"+name1,
                  defaultList
             );
    }
    else
    {
        if
        (
            (collisionType_ >= 2)
        and dictTransport.subDict("collisionData")
              .subDict("tabulatedInteractions")
              .subDict(collisionDataModel)
              .subDict("Omega22").found("shieldedCoulombPotential")
        )
        {
            isTabulatedData_ = false;
            
            if (collisionType_ == 3)
            {
                potentialType = "attractive";
            }
            else
            {
                potentialType = "repulsive";
            }
            
            cCDCoeffs2_[0] = dictTransport.subDict("collisionData")
                .subDict("tabulatedInteractions")
                .subDict(collisionDataModel).subDict("Omega22")
                .subDict("shieldedCoulombPotential").subDict(potentialType)
                .lookupOrDefault<scalar>("c2", 0.0);
            cCDCoeffs2_[1] = dictTransport.subDict("collisionData")
                .subDict("tabulatedInteractions")
                .subDict(collisionDataModel).subDict("Omega22")
                .subDict("shieldedCoulombPotential").subDict(potentialType)
                .lookupOrDefault<scalar>("C2", 0.0);
            cCDCoeffs2_[2] = dictTransport.subDict("collisionData")
                .subDict("tabulatedInteractions")
                .subDict(collisionDataModel).subDict("Omega22")
                .subDict("shieldedCoulombPotential").subDict(potentialType)
                .lookupOrDefault<scalar>("D2", 0.0);
        }
        else
        {
            FatalErrorIn
            (
                "collisionDataO::collisionDataO(const word&, const word&, "
                "const dictionary&, const dictionary&, const volScalarField&, "
                "const volScalarField&, const volScalarField&)"
            )   << collisionDataModel << " curve fit data missing for species "
                << "couple (" << name1 << ", " << name2 << ")."
                << exit(FatalError);
        }
    }
    
    if (isTabulatedData_)
    {
        piOmega1_[3] = exp(piOmega1_[3]);
        piOmega2_[3] = exp(piOmega2_[3]);
    }
}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::binaryDiffusionModels::collisionDataO::D() const
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
    const scalarField& Te = this->Te_.internalField();
    const scalarField& p = this->p_.internalField();
    const scalarField& pe = this->pe_.internalField();
    
    forAll(T, celli)
    {
        d[celli] = localkB_*T[celli]
            /(p[celli]*collisionTerm1(T[celli], pe[celli], Te[celli]));
    }
    
    forAll(this->T_.boundaryField(), patchi)
    {
        const fvPatchScalarField& pT = this->T_.boundaryField()[patchi];
        const fvPatchScalarField& pTe = this->Te_.boundaryField()[patchi];
	      const fvPatchScalarField& pp = this->p_.boundaryField()[patchi];
	      const fvPatchScalarField& ppe = this->pe_.boundaryField()[patchi];
	      
        fvPatchScalarField& pD = d.boundaryFieldRef()[patchi];

        forAll(pT, facei)
        {
            pD[facei] = localkB_*pT[facei]
                /(pp[facei]*collisionTerm1(pT[facei], ppe[facei], pTe[facei]));
        }
    }

    return tD;
}


Foam::tmp<Foam::scalarField> Foam::binaryDiffusionModels::collisionDataO::D
(
    const scalarField& p,
    const scalarField& pe,
    const scalarField& T,
    const scalarField& Te,
    const label patchi
) const
{
    tmp<scalarField> tD(new scalarField(T.size()));
    scalarField& d = tD.ref();

    forAll(T, facei)
    {
        d[facei] = localkB_*T[facei]
            /(p[facei]*collisionTerm1(T[facei], pe[facei], Te[facei]));
    }

    return tD;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
