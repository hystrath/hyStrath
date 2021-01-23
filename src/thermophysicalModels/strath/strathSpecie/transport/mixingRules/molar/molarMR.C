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

#include "molarMR.H"

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
    const volScalarField& T = thermo_.T();
    const volScalarField& p = thermo_.p();
    const volScalarField& CptrMix = thermo_.CptrMix();
    const volScalarField& CpvelMix = thermo_.CpvelMix();
    const volScalarField& CpMix = thermo_.CpMix();
    
    const scalarField& TCells = T.internalField();
    const scalarField& pCells = p.internalField();
    const scalarField& CpvelMixCells = CpvelMix.internalField();

    volScalarField& muMix = thermo_.mu();
    volScalarField& kappatrMix = thermo_.kappatr();
    volScalarField& kappaveMix = thermo_.kappave();
    volScalarField& kappaMix = thermo_.kappa();
    volScalarField& alphatrMix = thermo_.alphatr();
    volScalarField& alphaveMix = thermo_.alphave();
    volScalarField& alphaMix = thermo_.alpha();

    scalarField& muCells = muMix.primitiveFieldRef();
    scalarField& kappatrCells = kappatrMix.primitiveFieldRef();
    scalarField& kappaveCells = kappaveMix.primitiveFieldRef();
    scalarField& alphaveCells = alphaveMix.primitiveFieldRef();
    
    //- Initialisation
    muCells = 0.0;
    kappatrCells = 0.0;
    kappaveCells = 0.0;
        
    forAll(T.boundaryField(), patchi)
    {
        fvPatchScalarField& pmu = muMix.boundaryFieldRef()[patchi];
        fvPatchScalarField& pkappatr = kappatrMix.boundaryFieldRef()[patchi];
        fvPatchScalarField& pkappave = kappaveMix.boundaryFieldRef()[patchi];

        pmu = 0.0;
        pkappatr = 0.0;
        pkappave = 0.0;
    }

    forAll(species(), speciei)
    {
        const scalar R = thermo_.composition().R(speciei);
        
        scalar Rtr = R;
        scalar Rvel = 0.0;
        if (thermo_.composition().isElectron(speciei))
        {
            Rtr = 0.0;
            Rvel = R;
        }
        
        const volScalarField& Tve = thermo_.composition().Tv(speciei);
        const volScalarField& X = thermo_.composition().X(speciei);
        const volScalarField& Cvtr = thermo_.composition().Cvtr(speciei);
        const volScalarField& Cvvel = thermo_.composition().Cvvel(speciei);
        
        const scalarField& TveCells = Tve.internalField();
        const scalarField& XCells = X.internalField();
        const scalarField& CvtrCells = Cvtr.internalField();
        const scalarField& CvvelCells = Cvvel.internalField();

        scalarField& spmuCells = spmu_[speciei].primitiveFieldRef();
        scalarField& spkappatrCells = spkappatr_[speciei].primitiveFieldRef();
        scalarField& spkappaveCells = spkappave_[speciei].primitiveFieldRef();
        scalarField& spalphatrCells = spalphatr_[speciei].primitiveFieldRef();
        scalarField& spalphaveCells = spalphave_[speciei].primitiveFieldRef();

        //- Cell values
        forAll(XCells, celli)
        {
            spmuCells[celli] = mu(speciei, pCells[celli], TCells[celli]);
            muCells[celli] += XCells[celli]*spmuCells[celli];
            
            spkappatrCells[celli] =
                kappatr(speciei, pCells[celli], TCells[celli]);
            kappatrCells[celli] += XCells[celli]*spkappatrCells[celli];

            spkappaveCells[celli] =
                kappave(speciei, pCells[celli], TCells[celli], TveCells[celli]);
            kappaveCells[celli] += XCells[celli]*spkappaveCells[celli];
            
            spalphatrCells[celli] = spkappatrCells[celli]
               / max(CvtrCells[celli] + Rtr, 1e-12);
            spalphaveCells[celli] = spkappaveCells[celli]
               / max(CvvelCells[celli] + Rvel, 1e-12);
        }// end cells loop
        
        //- Patch values
        forAll(X.boundaryField(), patchi)
        {
            const fvPatchScalarField& pX = X.boundaryField()[patchi];
            const fvPatchScalarField& pp = p.boundaryField()[patchi];
            const fvPatchScalarField& pT = T.boundaryField()[patchi];
            const fvPatchScalarField& pTve = Tve.boundaryField()[patchi];
            const fvPatchScalarField& pCvtr = Cvtr.boundaryField()[patchi];
            const fvPatchScalarField& pCvvel = Cvvel.boundaryField()[patchi];

            fvPatchScalarField& pspmu =
                spmu_[speciei].boundaryFieldRef()[patchi];
            fvPatchScalarField& pmu = muMix.boundaryFieldRef()[patchi];
            fvPatchScalarField& pspkappatr =
                spkappatr_[speciei].boundaryFieldRef()[patchi];
            fvPatchScalarField& pkappatr =
                kappatrMix.boundaryFieldRef()[patchi];
            fvPatchScalarField& pspkappave =
                spkappave_[speciei].boundaryFieldRef()[patchi];
            fvPatchScalarField& pkappave =
                kappaveMix.boundaryFieldRef()[patchi];
            fvPatchScalarField& pspalphatr =
                spalphatr_[speciei].boundaryFieldRef()[patchi];
            fvPatchScalarField& pspalphave =
                spalphave_[speciei].boundaryFieldRef()[patchi];
                
            forAll(pX, facei)
            {
                pspmu[facei] = mu(speciei, pp[facei], pT[facei]);
                pmu[facei] += pX[facei]*pspmu[facei];
                
                pspkappatr[facei] = kappatr(speciei, pp[facei], pT[facei]);
                pkappatr[facei] += pX[facei]*pspkappatr[facei];
                
                pspkappave[facei] =
                    kappave(speciei, pp[facei], pT[facei], pTve[facei]);
                pkappave[facei] += pX[facei]*pspkappave[facei];
                
                pspalphatr[facei] = pspkappatr[facei]
                    /max(pCvtr[facei] + Rtr, 1e-12);
                pspalphave[facei] = pspkappave[facei]
                    /max(pCvvel[facei] + Rvel, 1e-12);
            }// end faces loop
        }// end patches loop
        
    }// end species loop

    kappaMix = kappatrMix + kappaveMix;
    alphatrMix = kappatrMix/CptrMix;
    alphaMix = kappaMix/CpMix;
    
    alphaveCells = kappaveCells/max(CpvelMixCells, 1e-12);

    forAll(T.boundaryField(), patchi)
    {
        const fvPatchScalarField& pCpvel = CpvelMix.boundaryField()[patchi];
        const fvPatchScalarField& pkappave =
            kappaveMix.boundaryFieldRef()[patchi];
            
        fvPatchScalarField& palphave = alphaveMix.boundaryFieldRef()[patchi];

        forAll(palphave, facei)
        {
            palphave[facei] = pkappave[facei]/max(pCpvel[facei], 1e-12);
        }
    }
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
            spalphatr_[speciei].write();
            spalphave_[speciei].write();
        }
    }

    if (writeKappaMixture_)
    {
        thermo_.kappatr().write();
        thermo_.kappave().write();
        thermo_.alpha().write();
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
