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

#include "WilkeMR.H"
#include "fvm.H"

//#include <time.h>

// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

template<class ThermoType>
void Foam::WilkeMR<ThermoType>::updatePhi()
{
    const volScalarField& p = thermo_.p();
    const volScalarField& Tt = thermo_.T();
    const scalarField& pCells = p.internalField();
    const scalarField& TtCells = Tt.internalField();

    forAll(species(), speciei)
    {
        volScalarField& phi = phi_[speciei];
        scalarField& phiCells = phi.primitiveFieldRef();

        phiCells = 0.0;

        forAll(Tt.boundaryField(), patchi)
        {
            phi.boundaryFieldRef()[patchi] = 0.0;
        }

        forAll(species(), speciej)
        {
            const volScalarField& Xj = thermo_.composition().X(speciej);
            const scalarField& XjCells = Xj.internalField();

            forAll(phiCells, celli)
            {
                if(miniXs_ < XjCells[celli])
                {
                    phiCells[celli] += XjCells[celli]*sqr
                    (
                        1.0 + sqrt(mu(speciei, pCells[celli], TtCells[celli])/mu(speciej, pCells[celli], TtCells[celli]))
                      * pow025(W(speciej)/W(speciei))
                    )
                      /sqrt(8.0*(1.0+W(speciei)/W(speciej)));
                }
            }

            forAll(Tt.boundaryField(), patchi)
            {
                const fvPatchScalarField& pXj = Xj.boundaryField()[patchi];
                const fvPatchScalarField& pp = p.boundaryField()[patchi];
                const fvPatchScalarField& pTt = Tt.boundaryField()[patchi];
                fvPatchScalarField& pphi = phi.boundaryFieldRef()[patchi];

                forAll(pTt, facei)
                {
                    if(miniXs_ < pXj[facei])
                    {
                        pphi[facei] += pXj[facei]*sqr
                        (
                            1.0 + sqrt(mu(speciei, pp[facei], pTt[facei])/mu(speciej, pp[facei], pTt[facei]))
                          * pow025(W(speciej)/W(speciei))
                        )
                         /sqrt(8.0*(1.0+W(speciei)/W(speciej)));
                    }
                }
            }
        }//end secondary species loop
    }//end primary species loop
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ThermoType>
Foam::WilkeMR<ThermoType>::WilkeMR
(
    rho2ReactionThermo& thermo,
    const compressibleTurbulenceModel& turbulence
)
:
    mixingRule(thermo, turbulence),

    miniXs_(1.0e-4)
{
    phi_.setSize(species().size());

    forAll(phi_, speciei)
    {
        phi_.set
        (
            speciei,
            new volScalarField
            (
                IOobject
                (
                    "phi_" + species()[speciei],
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh_,
                dimless
            )
        );
    }

    correct();
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class ThermoType>
void Foam::WilkeMR<ThermoType>::correct()
{
    /*std::clock_t start;
    std::clock_t restart;
    double duration;

    start = std::clock();*/
    updatePhi();
    /*duration = (std::clock() - start) / double(CLOCKS_PER_SEC);
    Info << "timer updatePhi() load: " << duration << endl;*/

    const volScalarField& p = thermo_.p();
    const volScalarField& Tt = thermo_.T();
    const scalarField& pCells = p.internalField();
    const scalarField& TtCells = Tt.internalField();

    volScalarField& muMix = thermo_.mu();
    volScalarField tempoMu = muMix;

    volScalarField& kappaMix = thermo_.kappatr();
    volScalarField& kappaveMix = thermo_.kappave();
    volScalarField tempoKappatr = kappaMix;
    volScalarField tempoKappave = kappaveMix;

    volScalarField& alphaMix = thermo_.alphatr();
    volScalarField& alphaveMix = thermo_.alphave();
    volScalarField tempoAlphatr = alphaMix;
    volScalarField tempoAlphave = alphaveMix;

    scalarField& muCells = tempoMu.primitiveFieldRef();
    scalarField& kappatrCells = tempoKappatr.primitiveFieldRef();
    scalarField& kappaveCells = tempoKappave.primitiveFieldRef();
    scalarField& alphatrCells = tempoAlphatr.primitiveFieldRef();
    scalarField& alphaveCells = tempoAlphave.primitiveFieldRef();

    //- Initialisations
    muCells = 0.0;
    kappatrCells = 0.0;
    kappaveCells = 0.0;
    alphatrCells = 0.0;
    alphaveCells = 0.0;

    forAll(tempoMu.boundaryField(), patchi)
    {
        fvPatchScalarField& pmu = tempoMu.boundaryFieldRef()[patchi];
        fvPatchScalarField& pkappatr = tempoKappatr.boundaryFieldRef()[patchi];
        fvPatchScalarField& pkappave = tempoKappave.boundaryFieldRef()[patchi];
        fvPatchScalarField& palphatr = tempoAlphatr.boundaryFieldRef()[patchi];
        fvPatchScalarField& palphave = tempoAlphave.boundaryFieldRef()[patchi];

        pmu = 0.0;
        pkappatr = 0.0;
        pkappave = 0.0;
        palphatr = 0.0;
        palphave = 0.0;

    }

    //restart = std::clock();
    //- Cell values
    forAll(species(), speciei)
    {
        const volScalarField& Tve = thermo_.composition().Tv(speciei);
        const volScalarField& X = thermo_.composition().X(speciei);

        const scalarField& TveCells = Tve.internalField();
        const scalarField& XCells = X.internalField();
        const scalarField& phiCells = phi_[speciei].internalField();
        scalarField& spmuCells = spmu_[speciei].primitiveFieldRef();
        scalarField& spkappatrCells = spkappatr_[speciei].primitiveFieldRef();
        scalarField& spkappaveCells = spkappave_[speciei].primitiveFieldRef();
        scalarField& spalphatrCells = spalphatr_[speciei].primitiveFieldRef();
        scalarField& spalphaveCells = spalphave_[speciei].primitiveFieldRef();

        forAll(XCells, celli)
        {
            const scalar factorCelli = XCells[celli]/phiCells[celli];

            spmuCells[celli] = mu(speciei, pCells[celli], TtCells[celli]);
            muCells[celli] += factorCelli*spmuCells[celli];

            spkappatrCells[celli] = kappatr(speciei, pCells[celli], TtCells[celli]);
            kappatrCells[celli] += factorCelli*spkappatrCells[celli];

            spkappaveCells[celli] = kappave(speciei, pCells[celli], TtCells[celli], TveCells[celli]);
            kappaveCells[celli] += factorCelli*spkappaveCells[celli];

            spalphatrCells[celli] = alphatr(speciei, pCells[celli], TtCells[celli]);
            alphatrCells[celli] += factorCelli*spalphatrCells[celli];

            spalphaveCells[celli] = alphave(speciei, pCells[celli], TtCells[celli], TveCells[celli]);
            alphaveCells[celli] += factorCelli*spalphaveCells[celli];
        }

        /*duration = (std::clock() - restart) / double(CLOCKS_PER_SEC);
        Info << "timer fields cell load: " << duration << endl;*/

        //- Patch values
        forAll(X.boundaryField(), patchi)
        {
            const fvPatchScalarField& pTve = Tve.boundaryField()[patchi];
            const fvPatchScalarField& pX = X.boundaryField()[patchi];

            const fvPatchScalarField& pTt = Tt.boundaryField()[patchi];
            const fvPatchScalarField& pp = p.boundaryField()[patchi];
            const fvPatchScalarField& pphi = phi_[speciei].boundaryField()[patchi];

            fvPatchScalarField& pspmu = spmu_[speciei].boundaryFieldRef()[patchi];
            fvPatchScalarField& pmu = tempoMu.boundaryFieldRef()[patchi];
            fvPatchScalarField& pspkappatr = spkappatr_[speciei].boundaryFieldRef()[patchi];
            fvPatchScalarField& pkappatr = tempoKappatr.boundaryFieldRef()[patchi];
            fvPatchScalarField& pspkappave = spkappave_[speciei].boundaryFieldRef()[patchi];
            fvPatchScalarField& pkappave = tempoKappave.boundaryFieldRef()[patchi];
            fvPatchScalarField& pspalphatr = spalphatr_[speciei].boundaryFieldRef()[patchi];
            fvPatchScalarField& palphatr = tempoAlphatr.boundaryFieldRef()[patchi];
            fvPatchScalarField& pspalphave = spalphave_[speciei].boundaryFieldRef()[patchi];
            fvPatchScalarField& palphave = tempoAlphave.boundaryFieldRef()[patchi];

            forAll(pX, facei)
            {
                const scalar pfactorFacei = pX[facei]/pphi[facei];

                pspmu[facei] = mu(speciei, pp[facei], pTt[facei]);
                pmu[facei] += pfactorFacei*pspmu[facei];

                pspkappatr[facei] = kappatr(speciei, pp[facei], pTt[facei]);
                pkappatr[facei] += pfactorFacei*pspkappatr[facei];

                pspkappave[facei] = kappave(speciei, pp[facei], pTt[facei], pTve[facei]);
                pkappave[facei] += pfactorFacei*pspkappave[facei];

                pspalphatr[facei] = alphatr(speciei, pp[facei], pTt[facei]);
                palphatr[facei] += pfactorFacei*pspalphatr[facei];

                pspalphave[facei] = alphave(speciei, pp[facei], pTt[facei], pTve[facei]);
                palphave[facei] += pfactorFacei*pspalphave[facei];
            }
        }

        /*duration = (std::clock() - restart) / double(CLOCKS_PER_SEC);
        Info << "timer fields bdry load: " << duration << endl;*/
    }

    muMix = tempoMu;
    kappaMix = tempoKappatr;
    kappaveMix = tempoKappave;
    alphaMix = tempoAlphatr;
    alphaveMix = tempoAlphave;
}


template<class ThermoType>
void Foam::WilkeMR<ThermoType>::write()
{
    if (writeMuSpecies_)
    {
        forAll(species(), speciei)
        {
            spmu_[speciei].write();
        }
    }

    if (writeMuMixture_)
    {
        thermo_.mu().write();
    }

    if (writeKappaSpecies_)
    {
        forAll(species(), speciei)
        {
            spkappatr_[speciei].write();
            spkappave_[speciei].write();
        }
    }

    if (writeKappaMixture_)
    {
        thermo_.kappatr().write();
        thermo_.kappave().write();
    }
}


template<class ThermoType>
bool Foam::WilkeMR<ThermoType>::read()
{
    if (regIOobject::read())
    {
        return true;
    }
    else
    {
        return false;
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
