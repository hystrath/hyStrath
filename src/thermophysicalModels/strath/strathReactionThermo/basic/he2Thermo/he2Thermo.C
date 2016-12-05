/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

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
    volScalarField::GeometricBoundaryField& hbf = h.boundaryField();

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
    const scalarField& TtCells = this->Tt_.internalField();
    scalarField& hetCells = this->het_.internalField();

    forAll(hetCells, celli)
    {
        hetCells[celli] =
            this->cellMixture(celli).HEt(pCells[celli], TtCells[celli]);
    }

    forAll(this->het_.boundaryField(), patchi)
    {
        this->het_.boundaryField()[patchi] = het
        (
            this->p_.boundaryField()[patchi],
            this->Tt_.boundaryField()[patchi],
            patchi
        );
    }
    
    this->hetBoundaryCorrection(this->het_);
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

// NEW VINCENT ****************************************************************
//-----                            het                                  -----//
template<class BasicThermo, class MixtureType>
Foam::tmp<Foam::volScalarField> Foam::he2Thermo<BasicThermo, MixtureType>::het
(
    const volScalarField& p,
    const volScalarField& Tt
) const
{
    const fvMesh& mesh = this->Tt_.mesh();

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
            this->het_.dimensions()
        )
    );

    volScalarField& het = thet();
    scalarField& hetCells = het.internalField();
    const scalarField& pCells = p.internalField();
    const scalarField& TtCells = Tt.internalField();

    forAll(hetCells, celli)
    {
        hetCells[celli] =
            this->cellMixture(celli).HEt(pCells[celli], TtCells[celli]);
    }

    forAll(het.boundaryField(), patchi)
    {
        scalarField& phet = het.boundaryField()[patchi];
        const scalarField& pp = p.boundaryField()[patchi];
        const scalarField& pTt = Tt.boundaryField()[patchi];

        forAll(phet, facei)
        {
            phet[facei] =
                this->patchFaceMixture(patchi, facei).HEt(pp[facei], pTt[facei]);
        }
    }

    return thet;
}


template<class BasicThermo, class MixtureType>
Foam::tmp<Foam::scalarField> Foam::he2Thermo<BasicThermo, MixtureType>::het
(
    const scalarField& p,
    const scalarField& Tt,
    const labelList& cells
) const
{
    tmp<scalarField> thet(new scalarField(Tt.size()));
    scalarField& het = thet();

    forAll(Tt, celli)
    {
        het[celli] = this->cellMixture(cells[celli]).HEt(p[celli], Tt[celli]);
    }

    return thet;
}


template<class BasicThermo, class MixtureType>
Foam::tmp<Foam::scalarField> Foam::he2Thermo<BasicThermo, MixtureType>::het
(
    const scalarField& p,
    const scalarField& Tt,
    const label patchi
) const
{
    tmp<scalarField> thet(new scalarField(Tt.size()));
    scalarField& het = thet();

    forAll(Tt, facei)
    {
        het[facei] =
            this->patchFaceMixture(patchi, facei).HEt(p[facei], Tt[facei]);
    }

    return thet;
}


//-----                            hev                                  -----//
template<class BasicThermo, class MixtureType>
Foam::tmp<Foam::volScalarField> Foam::he2Thermo<BasicThermo, MixtureType>::hev
(
    const volScalarField& p,
    const volScalarField& Tv
) const
{
    const fvMesh& mesh = this->Tt_.mesh();

    tmp<volScalarField> thev
    (
        new volScalarField
        (
            IOobject
            (
                "hev",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            this->hevMix_.dimensions()
        )
    );

    volScalarField& hev = thev();
    scalarField& hevCells = hev.internalField();
    const scalarField& pCells = p.internalField();
    const scalarField& TvCells = Tv.internalField();

    forAll(hevCells, celli)
    {
        hevCells[celli] =
            this->cellMixture(celli).HEv(pCells[celli], TvCells[celli]);
    }

    forAll(hev.boundaryField(), patchi)
    {
        scalarField& phev = hev.boundaryField()[patchi];
        const scalarField& pp = p.boundaryField()[patchi];
        const scalarField& pTv = Tv.boundaryField()[patchi];

        forAll(phev, facei)
        {
            phev[facei] =
                this->patchFaceMixture(patchi, facei).HEv(pp[facei], pTv[facei]);
        }
    }

    return thev;
}


template<class BasicThermo, class MixtureType>
Foam::tmp<Foam::scalarField> Foam::he2Thermo<BasicThermo, MixtureType>::hev
(
    const scalarField& p,
    const scalarField& Tv,
    const labelList& cells
) const
{
    tmp<scalarField> thev(new scalarField(Tv.size()));
    scalarField& hev = thev();

    forAll(Tv, celli)
    {
        hev[celli] = this->cellMixture(cells[celli]).HEv(p[celli], Tv[celli]);
    }

    return thev;
}


template<class BasicThermo, class MixtureType>
Foam::tmp<Foam::scalarField> Foam::he2Thermo<BasicThermo, MixtureType>::hev
(
    const scalarField& p,
    const scalarField& Tv,
    const label patchi
) const
{
    tmp<scalarField> thev(new scalarField(Tv.size()));
    scalarField& hev = thev();

    forAll(Tv, facei)
    {
        hev[facei] =
            this->patchFaceMixture(patchi, facei).HEv(p[facei], Tv[facei]);
    }

    return thev;
}


//-----                            heel                                 -----//
template<class BasicThermo, class MixtureType>
Foam::tmp<Foam::volScalarField> Foam::he2Thermo<BasicThermo, MixtureType>::heel
(
    const volScalarField& p,
    const volScalarField& Tv
) const
{
    const fvMesh& mesh = this->Tt_.mesh();

    tmp<volScalarField> theel
    (
        new volScalarField
        (
            IOobject
            (
                "heel",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            this->heelMix_.dimensions()
        )
    );

    volScalarField& heel = theel();
    scalarField& heelCells = heel.internalField();
    const scalarField& pCells = p.internalField();
    const scalarField& TvCells = Tv.internalField();

    forAll(heelCells, celli)
    {
        heelCells[celli] =
            this->cellMixture(celli).HEel(pCells[celli], TvCells[celli]);
    }

    forAll(heel.boundaryField(), patchi)
    {
        scalarField& pheel = heel.boundaryField()[patchi];
        const scalarField& pp = p.boundaryField()[patchi];
        const scalarField& pTv = Tv.boundaryField()[patchi];

        forAll(pheel, facei)
        {
            pheel[facei] =
                this->patchFaceMixture(patchi, facei).HEel(pp[facei], pTv[facei]);
        }
    }

    return theel;
}


template<class BasicThermo, class MixtureType>
Foam::tmp<Foam::scalarField> Foam::he2Thermo<BasicThermo, MixtureType>::heel
(
    const scalarField& p,
    const scalarField& Tv,
    const labelList& cells
) const
{
    tmp<scalarField> theel(new scalarField(Tv.size()));
    scalarField& heel = theel();

    forAll(Tv, celli)
    {
        heel[celli] = this->cellMixture(cells[celli]).HEel(p[celli], Tv[celli]);
    }

    return theel;
}


template<class BasicThermo, class MixtureType>
Foam::tmp<Foam::scalarField> Foam::he2Thermo<BasicThermo, MixtureType>::heel
(
    const scalarField& p,
    const scalarField& Tv,
    const label patchi
) const
{
    tmp<scalarField> theel(new scalarField(Tv.size()));
    scalarField& heel = theel();

    forAll(Tv, facei)
    {
        heel[facei] =
            this->patchFaceMixture(patchi, facei).HEel(p[facei], Tv[facei]);
    }

    return theel;
}


//-----                            hevel                                -----//
template<class BasicThermo, class MixtureType>
Foam::tmp<Foam::volScalarField> Foam::he2Thermo<BasicThermo, MixtureType>::hevel
(
    const volScalarField& p,
    const volScalarField& Tv
) const
{
    const fvMesh& mesh = this->Tt_.mesh();

    tmp<volScalarField> thevel
    (
        new volScalarField
        (
            IOobject
            (
                "hevel",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            hevel().dimensions()
        )
    );

    volScalarField& hevel = thevel();
    scalarField& hevelCells = hevel.internalField();
    const scalarField& pCells = p.internalField();
    const scalarField& TvCells = Tv.internalField();

    forAll(hevelCells, celli)
    {
        hevelCells[celli] =
            this->cellMixture(celli).HEvel(pCells[celli], TvCells[celli]);
    }

    forAll(hevel.boundaryField(), patchi)
    {
        scalarField& phevel = hevel.boundaryField()[patchi];
        const scalarField& pp = p.boundaryField()[patchi];
        const scalarField& pTv = Tv.boundaryField()[patchi];

        forAll(phevel, facei)
        {
            phevel[facei] =
                this->patchFaceMixture(patchi, facei).HEvel(pp[facei], pTv[facei]);
        }
    }

    return thevel;
}


template<class BasicThermo, class MixtureType>
Foam::tmp<Foam::scalarField> Foam::he2Thermo<BasicThermo, MixtureType>::hevel
(
    const scalarField& p,
    const scalarField& Tv,
    const labelList& cells
) const
{
    tmp<scalarField> thevel(new scalarField(Tv.size()));
    scalarField& hevel = thevel();

    forAll(Tv, celli)
    {
        hevel[celli] = this->cellMixture(cells[celli]).HEvel(p[celli], Tv[celli]);
    }

    return thevel;
}


template<class BasicThermo, class MixtureType>
Foam::tmp<Foam::scalarField> Foam::he2Thermo<BasicThermo, MixtureType>::hevel
(
    const scalarField& p,
    const scalarField& Tv,
    const label patchi
) const
{
    tmp<scalarField> thevel(new scalarField(Tv.size()));
    scalarField& hevel = thevel();

    forAll(Tv, facei)
    {
        hevel[facei] =
            this->patchFaceMixture(patchi, facei).HEvel(p[facei], Tv[facei]);
    }

    return thevel;
}


//-----                            zetar                                -----//
template<class BasicThermo, class MixtureType>
Foam::tmp<Foam::volScalarField> Foam::he2Thermo<BasicThermo, MixtureType>::zetar
(
    const volScalarField& p,
    const volScalarField& Tt,
    const volScalarField& Tv
) const
{
    const fvMesh& mesh = this->Tt_.mesh();

    tmp<volScalarField> tzetar
    (
        new volScalarField
        (
            IOobject
            (
                "zetar",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimless
        )
    );

    volScalarField& zetar = tzetar();
    scalarField& zetarCells = zetar.internalField();
    const scalarField& pCells = p.internalField();
    const scalarField& TtCells = Tt.internalField();
    const scalarField& TvCells = Tv.internalField();

    forAll(zetarCells, celli)
    {
        zetarCells[celli] =
            this->cellMixture(celli).zetar(pCells[celli], TtCells[celli], TvCells[celli]);
    }

    forAll(zetar.boundaryField(), patchi)
    {
        scalarField& pzetar = zetar.boundaryField()[patchi];
        const scalarField& pp = p.boundaryField()[patchi];
        const scalarField& pTt = Tt.boundaryField()[patchi];
        const scalarField& pTv = Tv.boundaryField()[patchi];

        forAll(pzetar, facei)
        {
            pzetar[facei] =
                this->patchFaceMixture(patchi, facei).zetar(pp[facei], pTt[facei], pTv[facei]);
        }
    }

    return tzetar;
}


template<class BasicThermo, class MixtureType>
Foam::tmp<Foam::scalarField> Foam::he2Thermo<BasicThermo, MixtureType>::zetar
(
    const scalarField& p,
    const scalarField& Tt,
    const scalarField& Tv,
    const labelList& cells
) const
{
    tmp<scalarField> tzetar(new scalarField(Tt.size()));
    scalarField& zetar = tzetar();

    forAll(zetar, celli)
    {
        zetar[celli] = this->cellMixture(cells[celli]).zetar(p[celli], Tt[celli], Tv[celli]);
    }

    return tzetar;
}


template<class BasicThermo, class MixtureType>
Foam::tmp<Foam::scalarField> Foam::he2Thermo<BasicThermo, MixtureType>::zetar
(
    const scalarField& p,
    const scalarField& Tt,
    const scalarField& Tv,
    const label patchi
) const
{
    tmp<scalarField> tzetar(new scalarField(Tt.size()));
    scalarField& zetar = tzetar();

    forAll(zetar, facei)
    {
        zetar[facei] =
            this->patchFaceMixture(patchi, facei).zetar(p[facei], Tt[facei], Tv[facei]);
    }

    return tzetar;
}


//-----                            zetav                                -----//
template<class BasicThermo, class MixtureType>
Foam::tmp<Foam::volScalarField> Foam::he2Thermo<BasicThermo, MixtureType>::zetav
(
    const volScalarField& p,
    const volScalarField& Tv
) const
{
    const fvMesh& mesh = this->Tt_.mesh();

    tmp<volScalarField> tzetav
    (
        new volScalarField
        (
            IOobject
            (
                "zetav",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimless
        )
    );

    volScalarField& zetav = tzetav();
    scalarField& zetavCells = zetav.internalField();
    const scalarField& pCells = p.internalField();
    const scalarField& TvCells = Tv.internalField();

    forAll(zetavCells, celli)
    {
        zetavCells[celli] =
            this->cellMixture(celli).zetav(pCells[celli], TvCells[celli]);
    }

    forAll(zetav.boundaryField(), patchi)
    {
        scalarField& pzetav = zetav.boundaryField()[patchi];
        const scalarField& pp = p.boundaryField()[patchi];
        const scalarField& pTv = Tv.boundaryField()[patchi];

        forAll(pzetav, facei)
        {
            pzetav[facei] =
                this->patchFaceMixture(patchi, facei).zetav(pp[facei], pTv[facei]);
        }
    }

    return tzetav;
}


template<class BasicThermo, class MixtureType>
Foam::tmp<Foam::scalarField> Foam::he2Thermo<BasicThermo, MixtureType>::zetav
(
    const scalarField& p,
    const scalarField& Tv,
    const labelList& cells
) const
{
    tmp<scalarField> tzetav(new scalarField(Tv.size()));
    scalarField& zetav = tzetav();

    forAll(zetav, celli)
    {
        zetav[celli] = this->cellMixture(cells[celli]).zetav(p[celli], Tv[celli]);
    }

    return tzetav;
}


template<class BasicThermo, class MixtureType>
Foam::tmp<Foam::scalarField> Foam::he2Thermo<BasicThermo, MixtureType>::zetav
(
    const scalarField& p,
    const scalarField& Tv,
    const label patchi
) const
{
    tmp<scalarField> tzetav(new scalarField(Tv.size()));
    scalarField& zetav = tzetav();

    forAll(zetav, facei)
    {
        zetav[facei] =
            this->patchFaceMixture(patchi, facei).zetav(p[facei], Tv[facei]);
    }

    return tzetav;
}


//-----                            zetael                                -----//
template<class BasicThermo, class MixtureType>
Foam::tmp<Foam::volScalarField> Foam::he2Thermo<BasicThermo, MixtureType>::zetael
(
    const volScalarField& p,
    const volScalarField& Tv
) const
{
    const fvMesh& mesh = this->Tt_.mesh();

    tmp<volScalarField> tzetael
    (
        new volScalarField
        (
            IOobject
            (
                "zetael",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimless
        )
    );

    volScalarField& zetael = tzetael();
    scalarField& zetaelCells = zetael.internalField();
    const scalarField& pCells = p.internalField();
    const scalarField& TvCells = Tv.internalField();

    forAll(zetaelCells, celli)
    {
        zetaelCells[celli] =
            this->cellMixture(celli).zetael(pCells[celli], TvCells[celli]);
    }

    forAll(zetael.boundaryField(), patchi)
    {
        scalarField& pzetael = zetael.boundaryField()[patchi];
        const scalarField& pp = p.boundaryField()[patchi];
        const scalarField& pTv = Tv.boundaryField()[patchi];

        forAll(pzetael, facei)
        {
            pzetael[facei] =
                this->patchFaceMixture(patchi, facei).zetael(pp[facei], pTv[facei]);
        }
    }

    return tzetael;
}


template<class BasicThermo, class MixtureType>
Foam::tmp<Foam::scalarField> Foam::he2Thermo<BasicThermo, MixtureType>::zetael
(
    const scalarField& p,
    const scalarField& Tv,
    const labelList& cells
) const
{
    tmp<scalarField> tzetael(new scalarField(Tv.size()));
    scalarField& zetael = tzetael();

    forAll(zetael, celli)
    {
        zetael[celli] = this->cellMixture(cells[celli]).zetael(p[celli], Tv[celli]);
    }

    return tzetael;
}


template<class BasicThermo, class MixtureType>
Foam::tmp<Foam::scalarField> Foam::he2Thermo<BasicThermo, MixtureType>::zetael
(
    const scalarField& p,
    const scalarField& Tv,
    const label patchi
) const
{
    tmp<scalarField> tzetael(new scalarField(Tv.size()));
    scalarField& zetael = tzetael();

    forAll(zetael, facei)
    {
        zetael[facei] =
            this->patchFaceMixture(patchi, facei).zetael(p[facei], Tv[facei]);
    }

    return tzetael;
}
// END NEW VINCENT ************************************************************


template<class BasicThermo, class MixtureType>
Foam::tmp<Foam::volScalarField>
Foam::he2Thermo<BasicThermo, MixtureType>::hc() const
{
    const fvMesh& mesh = this->Tt_.mesh();

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

    volScalarField& hcf = thc();
    scalarField& hcCells = hcf.internalField();

    forAll(hcCells, celli)
    {
        hcCells[celli] = this->cellMixture(celli).Hc();
    }

    forAll(hcf.boundaryField(), patchi)
    {
        scalarField& phc = hcf.boundaryField()[patchi];

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
    const scalarField& Tt,
    const scalarField& Tv,
    const label patchi
) const
{
    tmp<scalarField> tCp(new scalarField(Tt.size()));
    scalarField& cp = tCp();

    forAll(Tt, facei)
    {
        cp[facei] =
            this->patchFaceMixture(patchi, facei).Cp(p[facei], Tt[facei], Tv[facei]);
    }

    return tCp;
}


template<class BasicThermo, class MixtureType>
Foam::tmp<Foam::volScalarField>
Foam::he2Thermo<BasicThermo, MixtureType>::Cp() const
{
    const fvMesh& mesh = this->Tt_.mesh();

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

    volScalarField& cp = tCp();

    forAll(this->Tt_, celli)
    {
        cp[celli] =
            this->cellMixture(celli).Cp(this->p_[celli], this->Tt_[celli], this->Tv_[celli]);
    }

    forAll(this->Tt_.boundaryField(), patchi)
    {
        const fvPatchScalarField& pp = this->p_.boundaryField()[patchi];
        const fvPatchScalarField& pTt = this->Tt_.boundaryField()[patchi];
        const fvPatchScalarField& pTv = this->Tv_.boundaryField()[patchi];
        fvPatchScalarField& pCp = cp.boundaryField()[patchi];

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
    const scalarField& Tt,
    const scalarField& Tv,
    const label patchi
) const
{
    tmp<scalarField> tCv(new scalarField(Tt.size()));
    scalarField& cv = tCv();

    forAll(Tt, facei)
    {
        cv[facei] =
            this->patchFaceMixture(patchi, facei).Cv(p[facei], Tt[facei], Tv[facei]);
    }

    return tCv;
}


template<class BasicThermo, class MixtureType>
Foam::tmp<Foam::volScalarField>
Foam::he2Thermo<BasicThermo, MixtureType>::Cv() const
{
    const fvMesh& mesh = this->Tt_.mesh();

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

    volScalarField& cv = tCv();

    forAll(this->Tt_, celli)
    {
        cv[celli] =
            this->cellMixture(celli).Cv(this->p_[celli], this->Tt_[celli], this->Tv_[celli]);
    }

    forAll(this->T_.boundaryField(), patchi)
    {
        cv.boundaryField()[patchi] = Cv
        (
            this->p_.boundaryField()[patchi],
            this->Tt_.boundaryField()[patchi],
            this->Tv_.boundaryField()[patchi],
            patchi
        );
    }

    return tCv;
}


// BRAND NEW VINCENT **********************************************************
template<class BasicThermo, class MixtureType>
Foam::tmp<Foam::scalarField> Foam::he2Thermo<BasicThermo, MixtureType>::Cv_t
(
    const scalarField& p,
    const scalarField& Tt,
    const label patchi
) const
{
    tmp<scalarField> tCvt(new scalarField(Tt.size()));
    scalarField& Cvt = tCvt();

    forAll(Tt, facei)
    {
        Cvt[facei] =
            this->patchFaceMixture(patchi, facei).Cv_t(p[facei], Tt[facei]);
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

    volScalarField& Cvt = tCvt();

    forAll(this->Tt_, celli)
    {
       Cvt[celli] =
            //this->cellMixture_Cv_t(celli, this->p_[celli], this->Tt_[celli]);
            this->cellMixture(celli).Cv_t(this->p_[celli], this->Tt_[celli]);
    }

    forAll(this->Tt_.boundaryField(), patchi)
    {
        const fvPatchScalarField& pp = this->p_.boundaryField()[patchi];
        const fvPatchScalarField& pTt = this->Tt_.boundaryField()[patchi];
        fvPatchScalarField& pCvt = Cvt.boundaryField()[patchi];

        forAll(pTt, facei)
        {
            pCvt[facei] = /*this->patchFaceMixture_Cv_t
            (
                //&Foam::species::multiThermo<BasicThermo, MixtureType>::Cv_t,
                patchi, 
                facei,
                pp[facei],
                pTt[facei]
            );*/
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
    scalarField& Cvv = tCvv();

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
    const fvMesh& mesh = this->Tt_.mesh();

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

    volScalarField& Cvv = tCvv();

    forAll(this->Tv_, celli)
    {
       Cvv[celli] =
            this->cellMixture(celli).Cv_v(this->p_[celli], this->Tv_[celli]);
    }

    forAll(this->Tv_.boundaryField(), patchi)
    {
        const fvPatchScalarField& pp = this->p_.boundaryField()[patchi];
        const fvPatchScalarField& pTv = this->Tv_.boundaryField()[patchi];
        fvPatchScalarField& pCvv = Cvv.boundaryField()[patchi];

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
    const scalarField& Tt,
    const label patchi
) const
{
    tmp<scalarField> tCpt(new scalarField(Tt.size()));
    scalarField& Cpt = tCpt();

    forAll(Tt, facei)
    {
        Cpt[facei] =
            this->patchFaceMixture(patchi, facei).Cp_t(p[facei], Tt[facei]);
    }

    return tCpt;
}


template<class BasicThermo, class MixtureType>
Foam::tmp<Foam::volScalarField>
Foam::he2Thermo<BasicThermo, MixtureType>::Cp_t() const
{
    const fvMesh& mesh = this->Tt_.mesh();

    tmp<volScalarField> tCpt
    (
        new volScalarField
        (
            IOobject
            (
                "Cpt",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar("Cpt", dimEnergy/dimMass/dimTemperature, 0.0)
        )
    );

    volScalarField& Cpt = tCpt();

    forAll(this->Tt_, celli)
    {
       Cpt[celli] =
            //this->cellMixture_Cp_t(celli, this->p_[celli], this->Tt_[celli]);
            this->cellMixture(celli).Cp_t(this->p_[celli], this->Tt_[celli]);
    }

    forAll(this->Tt_.boundaryField(), patchi)
    {
        const fvPatchScalarField& pp = this->p_.boundaryField()[patchi];
        const fvPatchScalarField& pTt = this->Tt_.boundaryField()[patchi];
        fvPatchScalarField& pCpt = Cpt.boundaryField()[patchi];

        forAll(pTt, facei)
        {
            pCpt[facei] = /*this->patchFaceMixture_Cp_t
            (
                patchi, 
                facei,
                pp[facei],
                pTt[facei]
            );*/
            this->patchFaceMixture(patchi, facei).Cp_t
            (
                pp[facei],
                pTt[facei]
            );
        }
    }

    return tCpt;
}
// END BRAND NEW VINCENT ******************************************************


template<class BasicThermo, class MixtureType>
Foam::tmp<Foam::scalarField> Foam::he2Thermo<BasicThermo, MixtureType>::gamma
(
    const scalarField& p,
    const scalarField& Tt,
    const scalarField& Tv,
    const label patchi
) const
{
    tmp<scalarField> tgamma(new scalarField(Tt.size()));
    scalarField& gamma = tgamma();

    forAll(Tt, facei)
    {
        gamma[facei] =
            this->patchFaceMixture(patchi, facei).gamma(p[facei], Tt[facei], Tv[facei]);
    }

    return tgamma;
}


template<class BasicThermo, class MixtureType>
Foam::tmp<Foam::volScalarField>
Foam::he2Thermo<BasicThermo, MixtureType>::gamma() const
{
    const fvMesh& mesh = this->Tt_.mesh();

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

    volScalarField& gamma = tgamma();

    forAll(this->Tt_, celli)
    {
        gamma[celli] =
            this->cellMixture(celli).gamma(this->p_[celli], this->Tt_[celli], this->Tv_[celli]);
    }

    forAll(this->Tt_.boundaryField(), patchi)
    {
        const fvPatchScalarField& pp = this->p_.boundaryField()[patchi];
        const fvPatchScalarField& pTt = this->Tt_.boundaryField()[patchi];
        const fvPatchScalarField& pTv = this->Tv_.boundaryField()[patchi];
        fvPatchScalarField& pgamma = gamma.boundaryField()[patchi];

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


//NEW VINCENT *****************************************************************
template<class BasicThermo, class MixtureType>
Foam::tmp<Foam::scalarField> Foam::he2Thermo<BasicThermo, MixtureType>::TtHEt
(
    const scalarField& h,
    const scalarField& p,
    const scalarField& T0,
    const labelList& cells
) const
{
    tmp<scalarField> tTt(new scalarField(h.size()));
    scalarField& Tt = tTt();

    forAll(h, celli)
    {
        Tt[celli] =
            this->cellMixture(cells[celli]).TtHEt(h[celli], p[celli], T0[celli]);
    }

    return tTt;
}


template<class BasicThermo, class MixtureType>
Foam::tmp<Foam::scalarField> Foam::he2Thermo<BasicThermo, MixtureType>::TtHEt
(
    const scalarField& h,
    const scalarField& p,
    const scalarField& T0,
    const label patchi
) const
{

    tmp<scalarField> tTt(new scalarField(h.size()));
    scalarField& Tt = tTt();
    forAll(h, facei)
    {
        Tt[facei] = this->patchFaceMixture
        (
            patchi,
            facei
        ).TtHEt(h[facei], p[facei], T0[facei]);
    }

    return tTt;
}


template<class BasicThermo, class MixtureType>
Foam::tmp<Foam::scalarField> Foam::he2Thermo<BasicThermo, MixtureType>::TvHEv
(
    const scalarField& h,
    const scalarField& p,
    const scalarField& T0,
    const labelList& cells
) const
{
    tmp<scalarField> tTv(new scalarField(h.size()));
    scalarField& Tv = tTv();

    forAll(h, celli)
    {
        Tv[celli] =
            this->cellMixture(cells[celli]).TvHEv(h[celli], p[celli], T0[celli]);
    }

    return tTv;
}


template<class BasicThermo, class MixtureType>
Foam::tmp<Foam::scalarField> Foam::he2Thermo<BasicThermo, MixtureType>::TvHEv
(
    const scalarField& h,
    const scalarField& p,
    const scalarField& T0,
    const label patchi
) const
{

    tmp<scalarField> tTv(new scalarField(h.size()));
    scalarField& Tv = tTv();
    forAll(h, facei)
    {
        Tv[facei] = this->patchFaceMixture
        (
            patchi,
            facei
        ).TvHEv(h[facei], p[facei], T0[facei]);
    }

    return tTv;
}


template<class BasicThermo, class MixtureType>
Foam::tmp<Foam::scalarField> Foam::he2Thermo<BasicThermo, MixtureType>::TvelHEvel
(
    const scalarField& h,
    const scalarField& p,
    const scalarField& T0,
    const labelList& cells
) const
{
    tmp<scalarField> tTv(new scalarField(h.size()));
    scalarField& Tv = tTv();

    forAll(h, celli)
    {
        Tv[celli] =
            this->cellMixture(cells[celli]).TvelHEvel(h[celli], p[celli], T0[celli]);
    }

    return tTv;
}


template<class BasicThermo, class MixtureType>
Foam::tmp<Foam::scalarField> Foam::he2Thermo<BasicThermo, MixtureType>::TvelHEvel
(
    const scalarField& h,
    const scalarField& p,
    const scalarField& T0,
    const label patchi
) const
{

    tmp<scalarField> tTv(new scalarField(h.size()));
    scalarField& Tv = tTv();
    forAll(h, facei)
    {
        Tv[facei] = this->patchFaceMixture
        (
            patchi,
            facei
        ).TvelHEvel(h[facei], p[facei], T0[facei]);
    }

    return tTv;
}


template<class BasicThermo, class MixtureType>
Foam::tmp<Foam::scalarField> Foam::he2Thermo<BasicThermo, MixtureType>::alpha_t
(
    const scalarField& p,
    const scalarField& Tt,
    const label patchi
) const
{
    tmp<scalarField> talphat(new scalarField(Tt.size()));
    scalarField& alphat = talphat();

    forAll(Tt, facei)
    {
        alphat[facei] =
            this->patchFaceMixture(patchi, facei).alphatr(p[facei], Tt[facei]);
    }

    return talphat;
}


template<class BasicThermo, class MixtureType>
Foam::tmp<Foam::volScalarField>
Foam::he2Thermo<BasicThermo, MixtureType>::alpha_t() const
{
    const fvMesh& mesh = this->Tt_.mesh();

    tmp<volScalarField> talphat
    (
        new volScalarField
        (
            IOobject
            (
                "alpha_t",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar("alpha_t", dimMass/dimLength/dimTime, 0.0)
        )
    );

    volScalarField& alphat = talphat();

    forAll(this->Tt_, celli)
    {
       alphat[celli] =
            this->cellMixture(celli).alphatr(this->p_[celli], this->Tt_[celli]);
    }

    forAll(this->Tt_.boundaryField(), patchi)
    {
        const fvPatchScalarField& pp = this->p_.boundaryField()[patchi];
        const fvPatchScalarField& pTt = this->Tt_.boundaryField()[patchi];
        fvPatchScalarField& palphat = alphat.boundaryField()[patchi];

        forAll(pTt, facei)
        {
            palphat[facei] = this->patchFaceMixture(patchi, facei).alphatr(pp[facei], pTt[facei]);
        }
    }

    return talphat;
}


template<class BasicThermo, class MixtureType>
Foam::tmp<Foam::scalarField> Foam::he2Thermo<BasicThermo, MixtureType>::alpha_v
(
    const scalarField& p,
    const scalarField& Tv,
    const label patchi
) const
{
    tmp<scalarField> talphav(new scalarField(Tv.size()));
    scalarField& alphav = talphav();

    const fvPatchScalarField& Tt = this->Tt_.boundaryField()[patchi];
    forAll(Tv, facei)
    {
            alphav[facei] =
                this->patchFaceMixture(patchi, facei).alphave(p[facei], Tt[facei], Tv[facei]);
    }

    return talphav;
}


template<class BasicThermo, class MixtureType>
Foam::tmp<Foam::volScalarField>
Foam::he2Thermo<BasicThermo, MixtureType>::alpha_v() const
{
    const fvMesh& mesh = this->Tv_.mesh();

    tmp<volScalarField> talphav
    (
        new volScalarField
        (
            IOobject
            (
                "alpha_v",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar("alpha_v", dimMass/dimLength/dimTime, 0.0)
        )
    );

    volScalarField& alphav = talphav();

    forAll(this->Tv_, celli)
    {
        alphav[celli] =
            this->cellMixture(celli).alphave(this->p_[celli], this->Tt_[celli], this->Tv_[celli]);
    }

    forAll(this->Tv_.boundaryField(), patchi)
    {
        const fvPatchScalarField& pp = this->p_.boundaryField()[patchi];
        const fvPatchScalarField& pTt = this->Tt_.boundaryField()[patchi];
        const fvPatchScalarField& pTv = this->Tv_.boundaryField()[patchi];
        fvPatchScalarField& palphav = alphav.boundaryField()[patchi];

        forAll(pTv, facei)
        {
            palphav[facei] = this->patchFaceMixture(patchi, facei).alphave(pp[facei], pTt[facei], pTv[facei]);
        }
    }

    return talphav;
}
// END NEW VINCENT ******************************************************


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
