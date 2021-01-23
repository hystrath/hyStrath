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

#include "he2Thermo.H"

#include "gradient2TREnergyFvPatchScalarField.H"
#include "mixed2TREnergyFvPatchScalarField.H"

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class BasicThermo, class MixtureType>
void Foam::he2Thermo<BasicThermo, MixtureType>::
hetBoundaryCorrection(volScalarField& h)
{
    volScalarField::Boundary& hbf = h.boundaryFieldRef();

    forAll(hbf, patchi)
    {
        if (isA<gradient2TREnergyFvPatchScalarField>(hbf[patchi]))
        {
            refCast<gradient2TREnergyFvPatchScalarField>(hbf[patchi]).gradient()
                = hbf[patchi].fvPatchField::snGrad();
        }
        else if (isA<mixed2TREnergyFvPatchScalarField>(hbf[patchi]))
        {
            refCast<mixed2TREnergyFvPatchScalarField>(hbf[patchi]).refGrad()
                = hbf[patchi].fvPatchField::snGrad();
        }
    }
}


template<class BasicThermo, class MixtureType>
void Foam::he2Thermo<BasicThermo, MixtureType>::init()
{
    const scalarField& pCells = this->p_.internalField();
    const scalarField& TCells = this->T_.internalField();
    scalarField& hetCells = this->het().primitiveFieldRef();

    forAll(hetCells, celli)
    {
        hetCells[celli] =
            this->cellMixture(celli).HEt(pCells[celli], TCells[celli]);
    }

    forAll(this->het().boundaryField(), patchi)
    {
        this->het().boundaryFieldRef()[patchi] = het
        (
            this->p_.boundaryField()[patchi],
            this->T_.boundaryField()[patchi],
            patchi
        );
    }

    this->hetBoundaryCorrection(this->het());
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


template<class BasicThermo, class MixtureType>
Foam::he2Thermo<BasicThermo, MixtureType>::he2Thermo
(
    const fvMesh& mesh,
    const word& phaseName
)
:
    BasicThermo(mesh, phaseName),
    MixtureType(*this, mesh)
{
    init();
}


template<class BasicThermo, class MixtureType>
Foam::he2Thermo<BasicThermo, MixtureType>::he2Thermo
(
    const fvMesh& mesh,
    const dictionary& dict,
    const word& phaseName
)
:
    BasicThermo(mesh, dict, phaseName),
    MixtureType(*this, mesh)
{
    init();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class BasicThermo, class MixtureType>
Foam::he2Thermo<BasicThermo, MixtureType>::~he2Thermo()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicThermo, class MixtureType>
Foam::tmp<Foam::volScalarField> Foam::he2Thermo<BasicThermo, MixtureType>::het
(
    const volScalarField& p,
    const volScalarField& T
) const
{
    const fvMesh& mesh = this->T_.mesh();

    tmp<volScalarField> thet
    (
        new volScalarField
        (
            IOobject
            (
                "het",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            this->het().dimensions()
        )
    );

    volScalarField& het = thet.ref();
    scalarField& hetCells = het.primitiveFieldRef();
    const scalarField& pCells = p.internalField();
    const scalarField& TCells = T.internalField();

    forAll(hetCells, celli)
    {
        hetCells[celli] =
            this->cellMixture(celli).HEt(pCells[celli], TCells[celli]);
    }

    forAll(het.boundaryField(), patchi)
    {
        fvPatchScalarField& phet = het.boundaryFieldRef()[patchi];
        const fvPatchScalarField& pp = p.boundaryField()[patchi];
        const fvPatchScalarField& pT = T.boundaryField()[patchi];

        forAll(phet, facei)
        {
            phet[facei] =
                this->patchFaceMixture(patchi, facei).HEt(pp[facei], pT[facei]);
        }
    }

    return thet;
}


template<class BasicThermo, class MixtureType>
Foam::tmp<Foam::scalarField> Foam::he2Thermo<BasicThermo, MixtureType>::het
(
    const scalarField& p,
    const scalarField& T,
    const labelList& cells
) const
{
    tmp<scalarField> thet(new scalarField(T.size()));
    scalarField& het = thet.ref();

    forAll(T, celli)
    {
        het[celli] = this->cellMixture(cells[celli]).HEt(p[celli], T[celli]);
    }

    return thet;
}


template<class BasicThermo, class MixtureType>
Foam::tmp<Foam::scalarField> Foam::he2Thermo<BasicThermo, MixtureType>::het
(
    const scalarField& p,
    const scalarField& T,
    const label patchi
) const
{
    tmp<scalarField> thet(new scalarField(T.size()));
    scalarField& het = thet.ref();

    forAll(T, facei)
    {
        het[facei] =
            this->patchFaceMixture(patchi, facei).HEt(p[facei], T[facei]);
    }

    return thet;
}


template<class BasicThermo, class MixtureType>
Foam::tmp<Foam::volScalarField>
Foam::he2Thermo<BasicThermo, MixtureType>::hc() const
{
    const fvMesh& mesh = this->T_.mesh();

    tmp<volScalarField> thc
    (
        new volScalarField
        (
            IOobject
            (
                "hc",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimEnergy/dimMass
        )
    );

    volScalarField& hcf = thc.ref();
    scalarField& hcCells = hcf.primitiveFieldRef();

    forAll(hcCells, celli)
    {
        hcCells[celli] = this->cellMixture(celli).Hc();
    }

    forAll(hcf.boundaryField(), patchi)
    {
        fvPatchScalarField& phc = hcf.boundaryFieldRef()[patchi];

        forAll(phc, facei)
        {
            phc[facei] = this->patchFaceMixture(patchi, facei).Hc();
        }
    }

    return thc;
}


template<class BasicThermo, class MixtureType>
Foam::tmp<Foam::scalarField> Foam::he2Thermo<BasicThermo, MixtureType>::Cp
(
    const scalarField& p,
    const scalarField& T,
    const scalarField& Tv,
    const label patchi
) const
{
    tmp<scalarField> tCp(new scalarField(T.size()));
    scalarField& cp = tCp.ref();

    forAll(T, facei)
    {
        cp[facei] =
            this->patchFaceMixture(patchi, facei).Cp(p[facei], T[facei], Tv[facei]);
    }

    return tCp;
}


template<class BasicThermo, class MixtureType>
Foam::tmp<Foam::volScalarField>
Foam::he2Thermo<BasicThermo, MixtureType>::Cp() const
{
    const fvMesh& mesh = this->T_.mesh();

    tmp<volScalarField> tCp
    (
        new volScalarField
        (
            IOobject
            (
                "Cp",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimEnergy/dimMass/dimTemperature
        )
    );

    volScalarField& cp = tCp.ref();

    forAll(this->T_, celli)
    {
        cp[celli] =
            this->cellMixture(celli).Cp(this->p_[celli], this->T_[celli], this->Tv_[celli]);
    }

    forAll(this->T_.boundaryField(), patchi)
    {
        const fvPatchScalarField& pp = this->p_.boundaryField()[patchi];
        const fvPatchScalarField& pTt = this->T_.boundaryField()[patchi];
        const fvPatchScalarField& pTv = this->Tv_.boundaryField()[patchi];
        fvPatchScalarField& pCp = cp.boundaryFieldRef()[patchi];

        forAll(pTt, facei)
        {
            pCp[facei] =
                this->patchFaceMixture(patchi, facei).Cp(pp[facei], pTt[facei], pTv[facei]);
        }
    }

    return tCp;
}


template<class BasicThermo, class MixtureType>
Foam::tmp<Foam::scalarField>
Foam::he2Thermo<BasicThermo, MixtureType>::Cv
(
    const scalarField& p,
    const scalarField& T,
    const scalarField& Tv,
    const label patchi
) const
{
    tmp<scalarField> tCv(new scalarField(T.size()));
    scalarField& cv = tCv.ref();

    forAll(T, facei)
    {
        cv[facei] =
            this->patchFaceMixture(patchi, facei).Cv(p[facei], T[facei], Tv[facei]);
    }

    return tCv;
}


template<class BasicThermo, class MixtureType>
Foam::tmp<Foam::volScalarField>
Foam::he2Thermo<BasicThermo, MixtureType>::Cv() const
{
    const fvMesh& mesh = this->T_.mesh();

    tmp<volScalarField> tCv
    (
        new volScalarField
        (
            IOobject
            (
                "Cv",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimEnergy/dimMass/dimTemperature
        )
    );

    volScalarField& cv = tCv.ref();

    forAll(this->T_, celli)
    {
        cv[celli] =
            this->cellMixture(celli).Cv(this->p_[celli], this->T_[celli], this->Tv_[celli]);
    }

    forAll(this->T_.boundaryField(), patchi)
    {
        cv.boundaryFieldRef()[patchi] = Cv
        (
            this->p_.boundaryField()[patchi],
            this->T_.boundaryField()[patchi],
            this->Tv_.boundaryField()[patchi],
            patchi
        );
    }

    return tCv;
}


template<class BasicThermo, class MixtureType>
Foam::tmp<Foam::scalarField> Foam::he2Thermo<BasicThermo, MixtureType>::Cv_t
(
    const scalarField& p,
    const scalarField& T,
    const label patchi
) const
{
    tmp<scalarField> tCvt(new scalarField(T.size()));
    scalarField& Cvt = tCvt.ref();

    forAll(T, facei)
    {
        Cvt[facei] =
            this->patchFaceMixture(patchi, facei).Cv_t(p[facei], T[facei]);
    }

    return tCvt;
}


template<class BasicThermo, class MixtureType>
Foam::tmp<Foam::volScalarField>
Foam::he2Thermo<BasicThermo, MixtureType>::Cv_t() const
{
    const fvMesh& mesh = this->T_.mesh();

    tmp<volScalarField> tCvt
    (
        new volScalarField
        (
            IOobject
            (
                "Cvt",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar("Cvt", dimEnergy/dimMass/dimTemperature, 0.0)
        )
    );

    volScalarField& Cvt = tCvt.ref();

    forAll(this->T_, celli)
    {
       Cvt[celli] =
            this->cellMixture(celli).Cv_t(this->p_[celli], this->T_[celli]);
    }

    forAll(this->T_.boundaryField(), patchi)
    {
        const fvPatchScalarField& pp = this->p_.boundaryField()[patchi];
        const fvPatchScalarField& pTt = this->T_.boundaryField()[patchi];
        fvPatchScalarField& pCvt = Cvt.boundaryFieldRef()[patchi];

        forAll(pTt, facei)
        {
            pCvt[facei] = 
                this->patchFaceMixture(patchi, facei).Cv_t
                (
                    pp[facei],
                    pTt[facei]
                );
        }
    }

    return tCvt;
}


template<class BasicThermo, class MixtureType>
Foam::tmp<Foam::scalarField> Foam::he2Thermo<BasicThermo, MixtureType>::Cv_v
(
    const scalarField& p,
    const scalarField& Tv,
    const label patchi
) const
{
    tmp<scalarField> tCvv(new scalarField(Tv.size()));
    scalarField& Cvv = tCvv.ref();

    forAll(Tv, facei)
    {
        Cvv[facei] =
            this->patchFaceMixture(patchi, facei).Cv_v(p[facei], Tv[facei]);
    }

    return tCvv;
}


template<class BasicThermo, class MixtureType>
Foam::tmp<Foam::volScalarField>
Foam::he2Thermo<BasicThermo, MixtureType>::Cv_v() const
{
    const fvMesh& mesh = this->T_.mesh();

    tmp<volScalarField> tCvv
    (
        new volScalarField
        (
            IOobject
            (
                "Cvv",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar("Cvv", dimEnergy/dimMass/dimTemperature, 0.0)
        )
    );

    volScalarField& Cvv = tCvv.ref();

    forAll(this->Tv_, celli)
    {
       Cvv[celli] =
            this->cellMixture(celli).Cv_v(this->p_[celli], this->Tv_[celli]);
    }

    forAll(this->Tv_.boundaryField(), patchi)
    {
        const fvPatchScalarField& pp = this->p_.boundaryField()[patchi];
        const fvPatchScalarField& pTv = this->Tv_.boundaryField()[patchi];
        fvPatchScalarField& pCvv = Cvv.boundaryFieldRef()[patchi];

        forAll(pTv, facei)
        {
            pCvv[facei] = this->patchFaceMixture(patchi, facei).Cv_v
            (
                pp[facei],
                pTv[facei]
            );
        }
    }

    return tCvv;
}


template<class BasicThermo, class MixtureType>
Foam::tmp<Foam::scalarField> Foam::he2Thermo<BasicThermo, MixtureType>::Cp_t
(
    const scalarField& p,
    const scalarField& T,
    const label patchi
) const
{
    tmp<scalarField> tCpt(new scalarField(T.size()));
    scalarField& Cpt = tCpt.ref();

    forAll(T, facei)
    {
        Cpt[facei] =
            this->patchFaceMixture(patchi, facei).Cp_t(p[facei], T[facei]);
    }

    return tCpt;
}


//template<class BasicThermo, class MixtureType>
//Foam::tmp<Foam::volScalarField>
//Foam::he2Thermo<BasicThermo, MixtureType>::Cp_t() const
//{
//    const fvMesh& mesh = this->T_.mesh();

//    tmp<volScalarField> tCpt
//    (
//        new volScalarField
//        (
//            IOobject
//            (
//                "Cpt",
//                mesh.time().timeName(),
//                mesh,
//                IOobject::NO_READ,
//                IOobject::NO_WRITE
//            ),
//            mesh,
//            dimensionedScalar("Cpt", dimEnergy/dimMass/dimTemperature, 0.0)
//        )
//    );

//    volScalarField& Cpt = tCpt.ref();

//    forAll(this->T_, celli)
//    {
//       Cpt[celli] =
//            this->cellMixture(celli).Cp_t(this->p_[celli], this->T_[celli]);
//    }

//    forAll(this->T_.boundaryField(), patchi)
//    {
//        const fvPatchScalarField& pp = this->p_.boundaryField()[patchi];
//        const fvPatchScalarField& pTt = this->T_.boundaryField()[patchi];
//        fvPatchScalarField& pCpt = Cpt.boundaryFieldRef()[patchi];

//        forAll(pTt, facei)
//        {
//            pCpt[facei] =
//                this->patchFaceMixture(patchi, facei).Cp_t
//                (
//                    pp[facei],
//                    pTt[facei]
//                );
//        }
//    }

//    return tCpt;
//}


template<class BasicThermo, class MixtureType>
Foam::tmp<Foam::scalarField> Foam::he2Thermo<BasicThermo, MixtureType>::gamma
(
    const scalarField& p,
    const scalarField& T,
    const scalarField& Tv,
    const label patchi
) const
{
    tmp<scalarField> tgamma(new scalarField(T.size()));
    scalarField& gamma = tgamma.ref();

    forAll(T, facei)
    {
        gamma[facei] =
            this->patchFaceMixture(patchi, facei).gamma(p[facei], T[facei], Tv[facei]);
    }

    return tgamma;
}


template<class BasicThermo, class MixtureType>
Foam::tmp<Foam::volScalarField>
Foam::he2Thermo<BasicThermo, MixtureType>::gamma() const
{
    const fvMesh& mesh = this->T_.mesh();

    tmp<volScalarField> tgamma
    (
        new volScalarField
        (
            IOobject
            (
                "gamma",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimless
        )
    );

    volScalarField& gamma = tgamma.ref();

    forAll(this->T_, celli)
    {
        gamma[celli] =
            this->cellMixture(celli).gamma(this->p_[celli], this->T_[celli], this->Tv_[celli]);
    }

    forAll(this->T_.boundaryField(), patchi)
    {
        const fvPatchScalarField& pp = this->p_.boundaryField()[patchi];
        const fvPatchScalarField& pTt = this->T_.boundaryField()[patchi];
        const fvPatchScalarField& pTv = this->Tv_.boundaryField()[patchi];
        fvPatchScalarField& pgamma = gamma.boundaryFieldRef()[patchi];

        forAll(pTt, facei)
        {
            pgamma[facei] = this->patchFaceMixture(patchi, facei).gamma
            (
                pp[facei],
                pTt[facei],
                pTv[facei]
            );
        }
    }

    return tgamma;
}


//template<class BasicThermo, class MixtureType>
//Foam::tmp<Foam::scalarField> Foam::he2Thermo<BasicThermo, MixtureType>::alpha_t
//(
//    const scalarField& p,
//    const scalarField& T,
//    const label patchi
//) const
//{
//    tmp<scalarField> talphat(new scalarField(T.size()));
//    scalarField& alphat = talphat.ref();

//    forAll(T, facei)
//    {
//        alphat[facei] =
//            this->patchFaceMixture(patchi, facei).alphatr(p[facei], T[facei]);
//    }

//    return talphat;
//}


//template<class BasicThermo, class MixtureType>
//Foam::tmp<Foam::volScalarField>
//Foam::he2Thermo<BasicThermo, MixtureType>::alpha_t() const
//{
//    const fvMesh& mesh = this->T_.mesh();

//    tmp<volScalarField> talphat
//    (
//        new volScalarField
//        (
//            IOobject
//            (
//                "alpha_t",
//                mesh.time().timeName(),
//                mesh,
//                IOobject::NO_READ,
//                IOobject::NO_WRITE
//            ),
//            mesh,
//            dimensionedScalar("alpha_t", dimMass/dimLength/dimTime, 0.0)
//        )
//    );

//    volScalarField& alphat = talphat.ref();

//    forAll(this->T_, celli)
//    {
//       alphat[celli] =
//            this->cellMixture(celli).alphatr(this->p_[celli], this->T_[celli]);
//    }

//    forAll(this->T_.boundaryField(), patchi)
//    {
//        const fvPatchScalarField& pp = this->p_.boundaryField()[patchi];
//        const fvPatchScalarField& pT = this->T_.boundaryField()[patchi];
//        fvPatchScalarField& palphat = alphat.boundaryFieldRef()[patchi];

//        forAll(pT, facei)
//        {
//            palphat[facei] = this->patchFaceMixture(patchi, facei).alphatr(pp[facei], pT[facei]);
//        }
//    }

//    return talphat;
//}


//template<class BasicThermo, class MixtureType>
//Foam::tmp<Foam::scalarField> Foam::he2Thermo<BasicThermo, MixtureType>::alpha_v
//(
//    const scalarField& p,
//    const scalarField& Tv,
//    const label patchi
//) const
//{
//    tmp<scalarField> talphav(new scalarField(Tv.size()));
//    scalarField& alphav = talphav.ref();

//    const fvPatchScalarField& T = this->T_.boundaryField()[patchi];
//    forAll(Tv, facei)
//    {
//            alphav[facei] =
//                this->patchFaceMixture(patchi, facei).alphave(p[facei], T[facei], Tv[facei]);
//    }

//    return talphav;
//}


//template<class BasicThermo, class MixtureType>
//Foam::tmp<Foam::volScalarField>
//Foam::he2Thermo<BasicThermo, MixtureType>::alpha_v() const
//{
//    const fvMesh& mesh = this->Tv_.mesh();

//    tmp<volScalarField> talphav
//    (
//        new volScalarField
//        (
//            IOobject
//            (
//                "alpha_v",
//                mesh.time().timeName(),
//                mesh,
//                IOobject::NO_READ,
//                IOobject::NO_WRITE
//            ),
//            mesh,
//            dimensionedScalar("alpha_v", dimMass/dimLength/dimTime, 0.0)
//        )
//    );

//    volScalarField& alphav = talphav.ref();

//    forAll(this->Tv_, celli)
//    {
//        alphav[celli] =
//            this->cellMixture(celli).alphave
//            (
//                this->p_[celli],
//                this->T_[celli],
//                this->Tv_[celli]
//            );
//    }

//    forAll(this->Tv_.boundaryField(), patchi)
//    {
//        const fvPatchScalarField& pp = this->p_.boundaryField()[patchi];
//        const fvPatchScalarField& pT = this->T_.boundaryField()[patchi];
//        const fvPatchScalarField& pTv = this->Tv_.boundaryField()[patchi];
//        fvPatchScalarField& palphav = alphav.boundaryFieldRef()[patchi];

//        forAll(pTv, facei)
//        {
//            palphav[facei] =
//                this->patchFaceMixture(patchi, facei)
//                    .alphave(pp[facei], pT[facei], pTv[facei]);
//        }
//    }

//    return talphav;
//}


template<class BasicThermo, class MixtureType>
bool Foam::he2Thermo<BasicThermo, MixtureType>::read()
{
    if (BasicThermo::read())
    {
        MixtureType::read(*this);
        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
