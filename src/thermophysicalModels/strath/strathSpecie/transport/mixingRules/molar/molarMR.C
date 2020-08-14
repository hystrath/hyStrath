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

#include "molarMR.H"
#include "fvm.H"

// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ThermoType>
Foam::molarMR<ThermoType>::molarMR
(
    rho2ReactionThermo& thermo,
    const compressibleTurbulenceModel& turbulence
)
:
    mixingRule(thermo, turbulence),

    speciesThermo_
    (
        dynamic_cast<const multi2ComponentMixture<ThermoType>&>
            (this->thermo_).speciesData()
    )
{
    correct();
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class ThermoType>
void Foam::molarMR<ThermoType>::correct()
{
    const volScalarField& Tt = thermo_.T();
    const volScalarField& p = thermo_.p();

    const scalarField& TtCells = Tt.internalField();
    const scalarField& pCells = p.internalField();

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

    //- Cell values
    forAll(species(), speciei)
    {
        const volScalarField& Tve = thermo_.composition().Tv(speciei);
        const volScalarField& X = thermo_.composition().X(speciei);
        const scalarField& TveCells = Tve.internalField();
        const scalarField& XCells = X.internalField();

        scalarField& spmuCells = spmu_[speciei].primitiveFieldRef();
        scalarField& spkappatrCells = spkappatr_[speciei].primitiveFieldRef();
        scalarField& spkappaveCells = spkappave_[speciei].primitiveFieldRef();
        scalarField& spalphatrCells = spalphatr_[speciei].primitiveFieldRef();
        scalarField& spalphaveCells = spalphave_[speciei].primitiveFieldRef();

        forAll(XCells, celli)
        {
            spmuCells[celli] = mu(speciei, pCells[celli], TtCells[celli]);
            muCells[celli] += XCells[celli]*spmuCells[celli];
            spkappatrCells[celli] = kappatr(speciei, pCells[celli], TtCells[celli]);
            kappatrCells[celli] += XCells[celli]*spkappatrCells[celli];

            spkappaveCells[celli] = kappave(speciei, pCells[celli], TtCells[celli], TveCells[celli]);
            kappaveCells[celli] += XCells[celli]*spkappaveCells[celli];
            spalphatrCells[celli] = alphatr(speciei, pCells[celli], TtCells[celli]);
            alphatrCells[celli] += XCells[celli]*spalphatrCells[celli];
            spalphaveCells[celli] = alphave(speciei, pCells[celli], TtCells[celli], TveCells[celli]);
            alphaveCells[celli] += XCells[celli]*spalphaveCells[celli];
        }

        //- Patch values
        forAll(X.boundaryField(), patchi)
        {
            const fvPatchScalarField& pX = X.boundaryField()[patchi];
            const fvPatchScalarField& pp = p.boundaryField()[patchi];
            const fvPatchScalarField& pTt = Tt.boundaryField()[patchi];
            const fvPatchScalarField& pTve = Tve.boundaryField()[patchi];

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
                pspmu[facei] = mu(speciei, pp[facei], pTt[facei]);
                pmu[facei] += pX[facei]*pspmu[facei];
                pspkappatr[facei] = kappatr(speciei, pp[facei], pTt[facei]);
                pkappatr[facei] += pX[facei]*pspkappatr[facei];
                pspkappave[facei] = kappave(speciei, pp[facei], pTt[facei], pTve[facei]);
                pkappave[facei] += pX[facei]*pspkappave[facei];
                pspalphatr[facei] = alphatr(speciei, pp[facei], pTt[facei]);
                palphatr[facei] += pX[facei]*pspalphatr[facei];
                pspalphave[facei] = alphave(speciei, pp[facei], pTt[facei], pTve[facei]);
                palphave[facei] += pX[facei]*pspalphave[facei];
            }
        }
    }
    muMix = tempoMu;
    kappaMix = tempoKappatr;
    kappaveMix = tempoKappave;
    alphaMix = tempoAlphatr;
    alphaveMix = tempoAlphave;
}


template<class ThermoType>
void Foam::molarMR<ThermoType>::write()
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
            //spalphatr_[speciei].write();
            //spalphave_[speciei].write();
        }
    }

    if (writeKappaMixture_)
    {
        thermo_.kappatr().write();
        thermo_.kappave().write();
        //thermo_.alphatr().write();
        //thermo_.alphave().write();
    }
}


template<class ThermoType>
bool Foam::molarMR<ThermoType>::read()
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
