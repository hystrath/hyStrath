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

#include "AppletonBray.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ThermoType>
Foam::AppletonBray<ThermoType>::AppletonBray
(
    rho2ReactionThermo& thermo,
    const compressibleTurbulenceModel& turbulence
)
:
    relaxationTimeModelHE(thermo, turbulence),

    speciesThermo_
    (
        dynamic_cast<const multi2ComponentMixture<ThermoType>&>
            (this->thermo_).speciesData()
    ),

    electronListPosition_(species()["e-"]),
    sigma_er_(1.0e-20)
{
    const scalar ec = constant::electromagnetic::e.value();
    const scalar RR = constant::physicoChemical::R.value();
    const scalar NA = constant::physicoChemical::NA.value();
    const scalar pi = constant::mathematical::pi;
    const scalar kB = RR/NA;
    const scalar We = W(electronListPosition_);
    
    sigma_eIon_factor1_ = 8.0*pi*pow4(ec)/(27.0*sqr(kB));
    
    sigma_eIon_factor2_ = 9.0*pow3(kB)/(4.0*pi*pow6(ec));
    
    Qhe_factor_ = 3.0*RR*NA*sqrt(8.0*RR/(pi*We));
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class ThermoType>
void Foam::AppletonBray<ThermoType>::correct()
{
    const volScalarField& Tt = thermo_.Tt();
    const volScalarField& Te = thermo_.composition().Tv(electronListPosition_);
    const volScalarField& pDe = thermo_.composition().pD(electronListPosition_);
    const volScalarField& nDe = thermo_.composition().nD(electronListPosition_);

    const scalarField& TtCells = Tt.internalField();
    const scalarField& pDeCells = pDe.internalField();
    const scalarField& nDeCells = nDe.internalField();
    const scalarField& TeCells = Te.internalField();
    scalarField& QHECells = this->QHE_.primitiveFieldRef();

    QHECells = 0.0;

    forAll(Tt.boundaryField(), patchi)
    {
        fvPatchScalarField& pQHE = this->QHE_.boundaryFieldRef()[patchi];
        pQHE = 0.0;
    }

    forAll(species(), specier)
    {
        if (specier != electronListPosition_)
        {
            const volScalarField& pDr = thermo_.composition().pD(specier);
            const scalarField& pDrCells = pDr.internalField();

            if (speciesThermo_[specier].particleType() < 3)
            {
                forAll(pDrCells, celli)
                {
                    QHECells[celli] += sigma_er_*pDr[celli]/sqr(W(specier));
                }

                forAll(pDr.boundaryField(), patchi)
                {
                    const fvPatchScalarField& ppDr =
                        pDr.boundaryField()[patchi];
                    fvPatchScalarField& pQHE =
                        this->QHE_.boundaryFieldRef()[patchi];

                    forAll(ppDr, facei)
                    {
                        pQHE[facei] += sigma_er_*ppDr[facei]/sqr(W(specier));
                    }
                }
            }
            else
            {
                forAll(pDrCells, celli)
                {
                    scalar sigma_eIon = sigma_eIon_factor1_/sqr(TeCells[celli]);

                    if (nDeCells[celli] != 0.0)
                    {
                        sigma_eIon *=
                            log
                            (
                                1.0
                              + sigma_eIon_factor2_
                                  *pow3(TeCells[celli])/nDeCells[celli]
                            );
                    }

                    QHECells[celli] += sigma_eIon*pDr[celli]/sqr(W(specier));
                }

                forAll(pDr.boundaryField(), patchi)
                {
                    const fvPatchScalarField& pnDe =
                        nDe.boundaryField()[patchi];
                    const fvPatchScalarField& pTe = Te.boundaryField()[patchi];
                    const fvPatchScalarField& ppDr =
                        pDr.boundaryField()[patchi];
                    fvPatchScalarField& pQHE =
                        this->QHE_.boundaryFieldRef()[patchi];

                    forAll(ppDr, facei)
                    {
                        scalar sigma_eIon = sigma_eIon_factor1_/sqr(pTe[facei]);

                        if (pnDe[facei] != 0.0)
                        {
                            sigma_eIon *=
                                log
                                (
                                    1.0
                                 + sigma_eIon_factor2_
                                      *pow3(pTe[facei])/pnDe[facei]
                                );
                        }

                        pQHE[facei] += sigma_eIon*ppDr[facei]/sqr(W(specier));
                    }
                }
            }
        }
    }

    forAll(TtCells, celli)
    {
        Info<< "QHEC"<< tab << QHECells[celli]<<endl;
        Info<< "rhoe-"<< tab << pDeCells[celli]<<endl;
        QHECells[celli] *= Qhe_factor_*pDeCells[celli]*sqrt(TeCells[celli])
            *(TtCells[celli] - TeCells[celli]);
//        Info<< "QHEC"<< tab << QHECells[celli]<<endl;
    }

    forAll(Tt.boundaryField(), patchi)
    {
        const fvPatchScalarField& pTt = Tt.boundaryField()[patchi];
        const fvPatchScalarField& pTe = Te.boundaryField()[patchi];
        const fvPatchScalarField& ppDe = pDe.boundaryField()[patchi];
        fvPatchScalarField& pQHE = this->QHE_.boundaryFieldRef()[patchi];

        forAll(pTt, facei)
        {
            pQHE[facei] *= Qhe_factor_*ppDe[facei]*sqrt(pTe[facei])
                *(pTt[facei] - pTe[facei]);
        }
    }
}


template<class ThermoType>
bool Foam::AppletonBray<ThermoType>::read()
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
