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

#include "rho2ReactionThermo.H"
#include "fvMesh.H"

#include "wallFvPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(rho2ReactionThermo, 0);
    defineRunTimeSelectionTable(rho2ReactionThermo, fvMesh);
}

/* * * * * * * * * * * * * * * public static data * * * * * * * * * * * * * */

const Foam::scalar Foam::rho2ReactionThermo::vibrationalCutOffTemp = 150;

const Foam::scalar Foam::rho2ReactionThermo::miniYforSolvingEvEqn = 1.0e-4;

bool Foam::rho2ReactionThermo::temperatureFieldOutOfRange = false;

bool Foam::rho2ReactionThermo::hasCrashedButRecovered = false;

Foam::label Foam::rho2ReactionThermo::noCellsWithTemperatureFieldOutOfRange = 0;

Foam::FixedList<Foam::scalar, 2>
Foam::rho2ReactionThermo::minMaxTemperatureFieldOutOfRange;


// * * * * * * * * * * * * * private Member Functions  * * * * * * * * * * * //

Foam::word Foam::rho2ReactionThermo::transportToTypedef
(
    const word transportModel
)
{
    word typeDefName = word::null;

    if (transportModel == "constant")
    {
        typeDefName = "demConstGasEThermoPhysicsH2TGD";
    }
    else if (transportModel == "BlottnerEucken")
    {
        typeDefName = "demBEGasEThermoPhysicsH2TGD";
    }
    else if (transportModel == "powerLawEucken")
    {
        typeDefName = "demPLEGasEThermoPhysicsH2TGD";
    }
    else if (transportModel == "SutherlandEucken")
    {
        typeDefName = "demGasEThermoPhysicsH2TGD";
    }
    else if (transportModel == "CEA")
    {
        typeDefName = "demCEAGasEThermoPhysicsH2TGD";
    }

    return typeDefName;
}


void Foam::rho2ReactionThermo::correctChemFractions()
{
    const PtrList<Foam::volScalarField>& Y = composition().Y();
    const scalarField& pCells = this->p_.internalField();
    const scalarField& rhoCells = this->rho_.internalField();

    PtrList<Foam::volScalarField>& X = composition().X();
    PtrList<Foam::volScalarField>& nD = composition().nD();
    PtrList<Foam::volScalarField>& pP = composition().pP();
    PtrList<Foam::volScalarField>& pD = composition().pD();

    volScalarField& Wmix = composition().Wmix();

    scalarField sumX(0.0*X[0].internalField());
    scalarField sumnD(0.0*nD[0].internalField());
    scalarField sumpP(0.0*pP[0].internalField());
    scalarField sumpD(0.0*pD[0].internalField());

    forAll(Y, speciei)
    {
        const scalarField& YCells = Y[speciei].internalField();

        scalarField& XCells = X[speciei].primitiveFieldRef();
        scalarField& nDCells = nD[speciei].primitiveFieldRef();
        scalarField& pPCells = pP[speciei].primitiveFieldRef();
        scalarField& pDCells = pD[speciei].primitiveFieldRef();

        // This condition ensures that the sum of the chemical quantities are
        // bounded. Because Y is, quantities derived from Y will be.
        if (speciei < Y.size() - 1)
        {
            forAll(pCells, celli)
            {
                XCells[celli] =
                    composition().molarFraction(speciei, YCells[celli], celli);
                pPCells[celli] =
                    composition().partialPressure(XCells[celli], pCells[celli]);
                nDCells[celli] =
                    composition().numberDensity
                    (
                        speciei,
                        YCells[celli],
                        rhoCells[celli]
                    );
                pDCells[celli] =
                    composition().partialDensity
                    (
                        YCells[celli],
                        rhoCells[celli]
                    );

                sumX[celli] += XCells[celli];
                sumnD[celli] += nDCells[celli];
                sumpP[celli] += pPCells[celli];
                sumpD[celli] += pDCells[celli];
            }//end cells loop
        }
        else
        {
            forAll(pCells, celli)
            {
                XCells[celli] = max(1.0 - sumX[celli], 0.0);
                nDCells[celli] =
                    composition().numberDensity
                    (
                        speciei,
                        YCells[celli],
                        rhoCells[celli]
                    );
                pPCells[celli] = max(pCells[celli] - sumpP[celli], 0.0);
                pDCells[celli] = max(rhoCells[celli] - sumpD[celli], 0.0);
            }//end cells loop
        }
    }//end species loop


    forAll(this->p_.boundaryField(), patchi)
    {
        const fvPatchScalarField& pp  = this->p_.boundaryField()[patchi];
        const fvPatchScalarField& prho = this->rho_.boundaryField()[patchi];

        scalarField sumX(0.0*X[0].boundaryField()[patchi]);
        scalarField sumnD(0.0*nD[0].boundaryField()[patchi]);
        scalarField sumpP(0.0*pP[0].boundaryField()[patchi]);
        scalarField sumpD(0.0*pD[0].boundaryField()[patchi]);

        forAll(Y, speciei)
        {
            const fvPatchScalarField& pY = Y[speciei].boundaryField()[patchi];

            fvPatchScalarField& pX = X[speciei].boundaryFieldRef()[patchi];
            fvPatchScalarField& pnD = nD[speciei].boundaryFieldRef()[patchi];
            fvPatchScalarField& ppP = pP[speciei].boundaryFieldRef()[patchi];
            fvPatchScalarField& ppD = pD[speciei].boundaryFieldRef()[patchi];

            if (speciei < Y.size() - 1)
            {
                forAll(pp, facei)
                {
                    pX[facei] =
                        composition().molarFraction
                        (
                            speciei,
                            pY[facei],
                            patchi,
                            facei
                        );
                    pnD[facei] =
                        composition().numberDensity
                        (
                            speciei,
                            pY[facei],
                            prho[facei]
                        );
                    ppP[facei] =
                        composition().partialPressure(pX[facei], pp[facei]);
                    ppD[facei] =
                        composition().partialDensity(pY[facei], prho[facei]);

                    sumX[facei] += pX[facei];
                    sumnD[facei] += pnD[facei];
                    sumpP[facei] += ppP[facei];
                    sumpD[facei] += ppD[facei];
                }//end faces loop
            }
            else
            {
                forAll(pp, facei)
                {
                    pX[facei] = max(1.0 - sumX[facei], 0.0);
                    pnD[facei] =
                        composition().numberDensity
                        (
                            speciei,
                            pY[facei],
                            prho[facei]
                        );
                    ppP[facei] = max(pp[facei] - sumpP[facei], 0.0);
                    ppD[facei] = max(prho[facei] - sumpD[facei], 0.0);
                }//end faces loop
            }
        }//end species loop
    }//end patches loop
    
    // Calculation of the electron pressure 
    if (composition().contains("e-"))
    {
        this->pe_ = composition().pP("e-");
    }
    else
    {
        this->pe_.primitiveFieldRef() = 0.0;
        this->pe_.boundaryFieldRef() = 0.0;
    }
    
    Wmix = composition().molWeightMixture();
}


void Foam::rho2ReactionThermo::correctOverallTemperature()
{
    //- Declaration
    const PtrList<Foam::volScalarField>& Y = composition().Y();
    const PtrList<Foam::volScalarField>& X = composition().X();
    const PtrList<Foam::volScalarField>& Tv = composition().Tv();
    const scalarField& pCells = this->p_.internalField();
    const scalarField& TtCells = this->Tt_.internalField();

    PtrList<Foam::volScalarField>& zetar = composition().zetar();
    PtrList<Foam::volScalarField>& zetav = composition().zetav();
    PtrList<Foam::volScalarField>& zetael = composition().zetael();

    volScalarField numTMix = this->T_;
    volScalarField denTMix = this->T_;

    scalarField& zetarCellsMix = this->zetar_.primitiveFieldRef();
    scalarField& zetaelCellsMix = this->zetael_.primitiveFieldRef();
    scalarField& numTCellsMix = numTMix.primitiveFieldRef();
    scalarField& denTCellsMix = denTMix.primitiveFieldRef();
    scalarField& TCells = this->T_.primitiveFieldRef();

    //- Initialisation
    zetarCellsMix = 0.0;
    zetaelCellsMix = 0.0;
    numTCellsMix = 0.0;
    denTCellsMix = 0.0;

    forAll(this->p_.boundaryField(), patchi)
    {
        fvPatchScalarField& pzetarMix = this->zetar_.boundaryFieldRef()[patchi];
        fvPatchScalarField& pzetaelMix =
            this->zetael_.boundaryFieldRef()[patchi];
        fvPatchScalarField& pnumTMix = numTMix.boundaryFieldRef()[patchi];
        fvPatchScalarField& pdenTMix = denTMix.boundaryFieldRef()[patchi];

        pzetarMix = 0.0;
        pzetaelMix = 0.0;
        pnumTMix = 0.0;
        pdenTMix = 0.0;
    }

    //- Update degrees of freedom for each energy mode
    //- Field calculation
    forAll(Y, speciei)
    {
        //- Cells values
        const scalarField& YCells = Y[speciei].internalField();
        const scalarField& XCells = X[speciei].internalField();
        const scalarField& TvCells = Tv[speciei].internalField();
        const scalarField& zetavCells = zetav[speciei].internalField();
        scalarField& zetarCells = zetar[speciei].primitiveFieldRef();
        scalarField& zetaelCells = zetael[speciei].primitiveFieldRef();

        forAll(pCells, celli)
        {
            zetarCells[celli] =
                composition().zetar
                (
                    speciei,
                    pCells[celli],
                    TtCells[celli],
                    TvCells[celli]
                );
            zetaelCells[celli] =
                composition().zetael(speciei, pCells[celli], TvCells[celli]);
            zetarCellsMix[celli] += XCells[celli]*zetarCells[celli];
            zetaelCellsMix[celli] += XCells[celli]*zetaelCells[celli];

            if (composition().species()[speciei] == "e-")
            {
                numTCellsMix[celli] += 3.0*TvCells[celli]*YCells[celli];
                denTCellsMix[celli] += 3.0*YCells[celli];
            }
            else
            {
                numTCellsMix[celli] +=
                    (
                        (3.0 + zetarCells[celli])*TtCells[celli]
                      + (zetavCells[celli] + zetaelCells[celli])*TvCells[celli]
                    )*YCells[celli];
                denTCellsMix[celli] +=
                    (
                        3.0 + zetarCells[celli] + zetavCells[celli]
                      + zetaelCells[celli]
                    )*YCells[celli];
            }
        }//end cells loop

        //- Patch values calculation
        forAll(this->p_.boundaryField(), patchi)
        {
            const fvPatchScalarField& pp = this->p_.boundaryField()[patchi];
            const fvPatchScalarField& pTt = this->Tt_.boundaryField()[patchi];
            const fvPatchScalarField& pY = Y[speciei].boundaryField()[patchi];
            const fvPatchScalarField& pX = X[speciei].boundaryField()[patchi];
            const fvPatchScalarField& pTv = Tv[speciei].boundaryField()[patchi];
            const fvPatchScalarField& pzetav =
                zetav[speciei].boundaryField()[patchi];

            fvPatchScalarField& pzetar =
                zetar[speciei].boundaryFieldRef()[patchi];
            fvPatchScalarField& pzetael =
                zetael[speciei].boundaryFieldRef()[patchi];

            fvPatchScalarField& pzetarMix =
                this->zetar_.boundaryFieldRef()[patchi];
            fvPatchScalarField& pzetaelMix =
                this->zetael_.boundaryFieldRef()[patchi];
            fvPatchScalarField& pnumTMix = numTMix.boundaryFieldRef()[patchi];
            fvPatchScalarField& pdenTMix = denTMix.boundaryFieldRef()[patchi];

            forAll(pTt, facei)
            {
                pzetar[facei] =
                    composition().zetar
                    (
                        speciei,
                        pp[facei],
                        pTt[facei],
                        pTv[facei]
                    );
                pzetael[facei] =
                    composition().zetael(speciei, pp[facei], pTv[facei]);
                pzetarMix[facei] += pX[facei]*pzetar[facei];
                pzetaelMix[facei] += pX[facei]*pzetael[facei];

                if (composition().species()[speciei] == "e-")
                {
                    pnumTMix[facei] += 3.0*pTv[facei]*pY[facei];
                    pdenTMix[facei] += 3.0*pY[facei];
                }
                else
                {
                    pnumTMix[facei] +=
                        (
                            (3.0 + pzetar[facei])*pTt[facei] 
                          + (pzetav[facei] + pzetael[facei])*pTv[facei]
                        )*pY[facei];
                    pdenTMix[facei] +=
                        (
                            3.0 + pzetar[facei] + pzetav[facei] + pzetael[facei]
                        )*pY[facei];
                }
            }//end faces loop
        }//end patches loop
    }//end species loop

    //- Update overall temperature
    forAll(pCells, celli)
    {
        TCells[celli] = numTCellsMix[celli]/denTCellsMix[celli];
    }

    forAll(this->p_.boundaryField(), patchi)
    {
        const fvPatchScalarField& pnumTMix = numTMix.boundaryField()[patchi];
        const fvPatchScalarField& pdenTMix = denTMix.boundaryField()[patchi];
        fvPatchScalarField& pT = this->T_.boundaryFieldRef()[patchi];

        forAll(pT, facei)
        {
            pT[facei] = pnumTMix[facei]/pdenTMix[facei];
        }
    }
}


// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -//


void Foam::rho2ReactionThermo::initialiseLight()
{
    correctChemFractions();

    //- Declaration
    const bool downgradeSingleTemperature = this->downgradeSingleTemperature();
    const bool downgradeSingleTv = this->downgradeSingleTv();
    const PtrList<Foam::volScalarField>& Y = composition().Y();
    const PtrList<Foam::volScalarField>& X = composition().X();
    const scalarField& pCells = this->p_.internalField();
    const scalarField& TtCells = this->Tt_.internalField();

    PtrList<Foam::volScalarField>& hvel = composition().hevel();
    PtrList<Foam::volScalarField>& Tv = composition().Tv();

    const scalarField& TvCellsMix = this->Tv_.internalField();
    scalarField& htCellsMix = this->het_.primitiveFieldRef();
    scalarField& hvelCellsMix = this->hevel().primitiveFieldRef();
    scalarField& hCellsMix = composition().e().primitiveFieldRef();
    scalarField& psiCellsMix = this->psi_.primitiveFieldRef();
    scalarField& rhoCellsMix = this->rho_.primitiveFieldRef();

    //- Cells values
    forAll(pCells, celli)
    {
        //- Initialisation
        htCellsMix[celli] = 0.0;
        hvelCellsMix[celli] = 0.0;
        hCellsMix[celli] = 0.0;
        psiCellsMix[celli] = 0.0;
        rhoCellsMix[celli] = 0.0;

        //- Calculation
        forAll(Y, speciei)
        {
            const scalarField& YCells = Y[speciei].internalField();
            htCellsMix[celli] += YCells[celli]*composition().HEt(speciei, pCells[celli], TtCells[celli]);
        }

        forAll(Y, speciei)
        {
            const scalarField& YCells = Y[speciei].internalField();
            const scalarField& XCells = X[speciei].internalField();

            scalarField& hvelCells = hvel[speciei].primitiveFieldRef();
            scalarField& TvCells = Tv[speciei].primitiveFieldRef();

            if (TvCells[celli] != 0.0)
            {
                if (downgradeSingleTemperature)
                {
                    TvCells[celli] = TtCells[celli];
                }
                else if (downgradeSingleTv)
                {
                    TvCells[celli] = TvCellsMix[celli];
                }

                hvelCells[celli] = composition().HEvel(speciei, pCells[celli], TvCells[celli]);
                hvelCellsMix[celli] += YCells[celli]*hvelCells[celli];

                hCellsMix[celli] += YCells[celli]*hvelCells[celli];
            }

            if (composition().particleType(speciei) > 0)
            {
                // NOTE VINCENT: Because of the way composition().rho is implemented, it requires X instead of Y
                psiCellsMix[celli] += XCells[celli]*composition().psi(speciei, pCells[celli], TtCells[celli]);
                rhoCellsMix[celli] += XCells[celli]*composition().rho(speciei, pCells[celli], TtCells[celli]);
            }
            else
            {
                psiCellsMix[celli] += XCells[celli]*composition().psi(speciei, pCells[celli], TvCells[celli]);
                rhoCellsMix[celli] += XCells[celli]*composition().rho(speciei, pCells[celli], TvCells[celli]);
            }
        }

        hCellsMix[celli] += htCellsMix[celli];
    }

    //- Patch values calculations
    forAll(this->Tt_.boundaryField(), patchi)
    {
        const fvPatchScalarField& pp = this->p_.boundaryField()[patchi];
        const fvPatchScalarField& pTt = this->Tt_.boundaryField()[patchi];
        fvPatchScalarField& pTvMix = this->Tv_.boundaryFieldRef()[patchi];

        fvPatchScalarField& phtMix = this->het_.boundaryFieldRef()[patchi];
        fvPatchScalarField& phvelMix = this->hevel().boundaryFieldRef()[patchi];
        fvPatchScalarField& phMix = composition().e().boundaryFieldRef()[patchi];

        if (pTt.fixesValue())
        // the temperature is fixed ... the energy is calculated at patches
        {
            forAll(pTt, facei)
            {
                //- Initialisation
                phtMix[facei] = 0.0;
                phvelMix[facei] = 0.0;
                phMix[facei] = 0.0;

                //- Calculation
                forAll(Y, speciei)
                {
                    const fvPatchScalarField& pY = Y[speciei].boundaryField()[patchi];
                    fvPatchScalarField& pTv = Tv[speciei].boundaryFieldRef()[patchi];

                    fvPatchScalarField& phvel = hvel[speciei].boundaryFieldRef()[patchi];

                    phtMix[facei] += pY[facei]*composition().HEt(speciei, pp[facei], pTt[facei]);

                    if (pTv[facei] != 0.0)
                    {
                        if (downgradeSingleTv)
                        {
                            pTv[facei] = pTvMix[facei];
                        }

                        phvel[facei] = composition().HEvel(speciei, pp[facei], pTv[facei]);
                        phvelMix[facei] += pY[facei]*phvel[facei];

                        phMix[facei] += pY[facei]*phvel[facei];
                    }
                }

                phMix[facei] += phtMix[facei];
            }
        }
        else
        // condition on the energy fields ... the temperatures are calculated at patches
        {
            forAll(pTt, facei)
            {
                //- Initialisation
                phtMix[facei] = 0.0;
                phvelMix[facei] = 0.0; // NEW VINCENT 14/02/2017
                phMix[facei] = 0.0;

                //- Calculation
                forAll(Y, speciei)
                {
                    const fvPatchScalarField& pY = Y[speciei].boundaryField()[patchi];
                    fvPatchScalarField& phvel = hvel[speciei].boundaryFieldRef()[patchi];
                    fvPatchScalarField& pTv = Tv[speciei].boundaryFieldRef()[patchi];

                    phtMix[facei] += pY[facei]*composition().HEt(speciei, pp[facei], pTt[facei]);

                    if (pTv[facei] != 0.0)
                    {
                        if (downgradeSingleTemperature)
                        {
                            phvel[facei] = composition().HEvel(speciei, pp[facei], pTt[facei]); // NEW VINCENT 14/02/2017
                            pTv[facei] = pTt[facei];
                        }
                        else if (downgradeSingleTv)
                        {
                            phvel[facei] = composition().HEvel(speciei, pp[facei], pTvMix[facei]); // NEW VINCENT 14/02/2017
                            pTv[facei] = pTvMix[facei];
                        }
                        else if (composition().vibTempAssociativity(speciei) == -1)
                        {
                            phvel[facei] = composition().HEvel(speciei, pp[facei], pTv[facei]); // NEW VINCENT 14/02/2017
                            pTv[facei] = TvelEvels(speciei, phvel[facei], pp[facei], pTv[facei]);
                            //composition().TvelHEvel(speciei, phvel[facei], pp[facei], pTv[facei]);
                        }
                        else
                        {
                            const fvPatchScalarField& pTvMol = Tv[composition().vibTempAssociativity(speciei)].boundaryField()[patchi];
                            phvel[facei] = composition().HEvel(speciei, pp[facei], pTvMol[facei]); // NEW VINCENT 14/02/2017
                            pTv[facei] = pTvMol[facei];
                        }

                        phvelMix[facei] += pY[facei]*phvel[facei];

                        phMix[facei] += pY[facei]*phvel[facei];
                    }
                }//end species loop

                phMix[facei] += phtMix[facei];
            }//end faces loop
        }

        fvPatchScalarField& ppsiMix = this->psi_.boundaryFieldRef()[patchi];
        fvPatchScalarField& prhoMix = this->rho_.boundaryFieldRef()[patchi];

        //- Initialisation
        forAll(pTt, facei)
        {
            ppsiMix[facei] = 0.0;
            prhoMix[facei] = 0.0;
        }

        //- Calculation
        forAll(Y, speciei)
        {
            const fvPatchScalarField& pX = X[speciei].boundaryField()[patchi];
            const fvPatchScalarField& pTv = Tv[speciei].boundaryField()[patchi];

            forAll(pTt, facei)
            {
                if (composition().particleType(speciei) > 0)
                {
                    ppsiMix[facei] += pX[facei]*composition().psi(speciei, pp[facei], pTt[facei]);
                    prhoMix[facei] += pX[facei]*composition().rho(speciei, pp[facei], pTt[facei]);
                }
                else
                {
                    ppsiMix[facei] += pX[facei]*composition().psi(speciei, pp[facei], pTv[facei]);
                    prhoMix[facei] += pX[facei]*composition().rho(speciei, pp[facei], pTv[facei]);
                }
            }//end faces loop
        }//end species loop
    }
}


// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -//


void Foam::rho2ReactionThermo::initialise()
{
    correctChemFractions();

    //- Declarations
    const bool downgradeSingleTemperature = this->downgradeSingleTemperature();
    const bool downgradeSingleTv = this->downgradeSingleTv();
    const PtrList<Foam::volScalarField>& Y = composition().Y();
    const PtrList<Foam::volScalarField>& X = composition().X();
    const scalarField& pCells = this->p_.internalField();
    const scalarField& TtCells = this->Tt_.internalField();

    PtrList<Foam::volScalarField>& hv = composition().hev();
    PtrList<Foam::volScalarField>& hel = composition().heel();
    PtrList<Foam::volScalarField>& hvel = composition().hevel();
    PtrList<Foam::volScalarField>& Tv = composition().Tv();
    PtrList<Foam::volScalarField>& zetav = composition().zetav();

    //PtrList<PtrList<Foam::volScalarField> >& hvel_mode = composition().hevel_mode(); // NEW VINCENT 23/03/2016 TODO ONGOING WORK
    //PtrList<PtrList<Foam::volScalarField> >& Tv_mode = composition().Tv_mode(); // NEW VINCENT 14/03/2016 TODO ONGOING WORK
    //PtrList<PtrList<Foam::volScalarField> >& zetav_mode = composition().zetav_mode(); // NEW VINCENT 14/03/2016 TODO ONGOING WORK

    scalarField& TvCellsMix = this->Tv_.primitiveFieldRef();
    scalarField& htCellsMix = this->het_.primitiveFieldRef();
    scalarField& hvCellsMix = this->hevMix_.primitiveFieldRef();
    scalarField& helCellsMix = this->heelMix_.primitiveFieldRef();
    scalarField& hvelCellsMix = this->hevel().primitiveFieldRef();
    scalarField& hCellsMix = composition().e().primitiveFieldRef();
    scalarField& zetavCellsMix = this->zetav_.primitiveFieldRef();
    scalarField totZetavCellsMix = zetavCellsMix;
    scalarField& psiCellsMix = this->psi_.primitiveFieldRef();
    scalarField& rhoCellsMix = this->rho_.primitiveFieldRef();

    //- Cells values
    forAll(pCells, celli)
    {
        //- Initialisation
        htCellsMix[celli] = 0.0;
        hvCellsMix[celli] = 0.0;
        helCellsMix[celli] = 0.0;
        hvelCellsMix[celli] = 0.0;
        hCellsMix[celli] = 0.0;
        if (not downgradeSingleTv)
        {
            TvCellsMix[celli] = 0.0;
        }
        zetavCellsMix[celli] = 0.0;
        totZetavCellsMix[celli] = 0.0;

        psiCellsMix[celli] = 0.0;
        rhoCellsMix[celli] = 0.0;

        //- Calculation
        forAll(Y, speciei)
        {
            const scalarField& YCells = Y[speciei].internalField();
            htCellsMix[celli] += YCells[celli]*composition().HEt(speciei, pCells[celli], TtCells[celli]);
        }

        forAll(Y, speciei)
        {
            const scalarField& YCells = Y[speciei].internalField();
            const scalarField& XCells = X[speciei].internalField();

            scalarField& hvCells = hv[speciei].primitiveFieldRef();
            scalarField& helCells = hel[speciei].primitiveFieldRef();
            scalarField& hvelCells = hvel[speciei].primitiveFieldRef();
            scalarField& TvCells = Tv[speciei].primitiveFieldRef();
            scalarField& zetavCells = zetav[speciei].primitiveFieldRef();

            if (TvCells[celli] != 0.0)
            {
                if (downgradeSingleTemperature)
                {
                    TvCells[celli] = TtCells[celli];
                }
                else if (downgradeSingleTv)
                {
                    TvCells[celli] = TvCellsMix[celli];
                }

                hvCells[celli] = composition().HEv(speciei, pCells[celli], TvCells[celli]);
                helCells[celli] = composition().HEel(speciei, pCells[celli], TvCells[celli]);
                hvelCells[celli] = composition().HEvel(speciei, pCells[celli], TvCells[celli]);
                /*forAll(hvel_mode[speciei], vibMode) TODO ONGOING WORK
                {
                    hvel_mode[speciei][vibMode].internalField()[celli] =
                        composition().HEvel_mode
                        (
                            speciei,
                            vibMode,
                            pCells[celli],
                            Tv_mode[speciei][vibMode].internalField()[celli]
                        ); // NEW VINCENT 23/03/2016
                }*/

                zetavCells[celli] = composition().zetav(speciei, pCells[celli], TvCells[celli]);

                hvCellsMix[celli] += YCells[celli]*hvCells[celli];
                helCellsMix[celli] += YCells[celli]*helCells[celli];
                hvelCellsMix[celli] += YCells[celli]*hvelCells[celli];
                hCellsMix[celli] += YCells[celli]*hvelCells[celli];
                if (not downgradeSingleTv)
                {
                    TvCellsMix[celli] += XCells[celli]*zetavCells[celli]*TvCells[celli];
                }
                zetavCellsMix[celli] += XCells[celli]*zetavCells[celli];
                totZetavCellsMix[celli] += XCells[celli]*zetavCells[celli];
            }

            if (composition().particleType(speciei) > 0)
            {
                // NOTE VINCENT: Because of the way composition().rho is implemented, it requires X instead of Y
                psiCellsMix[celli] += XCells[celli]*composition().psi(speciei, pCells[celli], TtCells[celli]);
                rhoCellsMix[celli] += XCells[celli]*composition().rho(speciei, pCells[celli], TtCells[celli]);
            }
            else
            {
                psiCellsMix[celli] += XCells[celli]*composition().psi(speciei, pCells[celli], TvCells[celli]);
                rhoCellsMix[celli] += XCells[celli]*composition().rho(speciei, pCells[celli], TvCells[celli]);
            }
        }//end species loop

        if (not downgradeSingleTv and totZetavCellsMix[celli] != 0.0)
        {
            TvCellsMix[celli] /= totZetavCellsMix[celli];
        }

        hCellsMix[celli] += htCellsMix[celli];
    }//end cells loop

    //- Patch values
    forAll(this->Tt_.boundaryField(), patchi)
    {
        const fvPatchScalarField& pp = this->p_.boundaryField()[patchi];
        const fvPatchScalarField& pTt = this->Tt_.boundaryField()[patchi];

        fvPatchScalarField& phtMix = this->het_.boundaryFieldRef()[patchi];
        fvPatchScalarField& phvMix = this->hevMix_.boundaryFieldRef()[patchi];
        fvPatchScalarField& phelMix = this->heelMix_.boundaryFieldRef()[patchi];
        fvPatchScalarField& phvelMix = this->hevel().boundaryFieldRef()[patchi];
        fvPatchScalarField& phMix = composition().e().boundaryFieldRef()[patchi];
        fvPatchScalarField& pTvMix = this->Tv_.boundaryFieldRef()[patchi];
        fvPatchScalarField& pzetavMix = this->zetav_.boundaryFieldRef()[patchi];
        fvPatchScalarField ptotZetavMix = pzetavMix;

        if (pTt.fixesValue())
        // the temperature is fixed ... the energy is calculated at patches
        {
            forAll(pTt, facei)
            {
                //- Initialisation
                phtMix[facei] = 0.0;
                phvMix[facei] = 0.0;
                phelMix[facei] = 0.0;
                phvelMix[facei] = 0.0;
                phMix[facei] = 0.0;
                if (not downgradeSingleTv)
                {
                    pTvMix[facei] = 0.0;
                }
                pzetavMix[facei] = 0.0;
                ptotZetavMix[facei] = 0.0;

                //- Calculation
                forAll(Y, speciei)
                {
                    const fvPatchScalarField& pY = Y[speciei].boundaryField()[patchi];
                    const fvPatchScalarField& pX = X[speciei].boundaryField()[patchi];

                    fvPatchScalarField& phv = hv[speciei].boundaryFieldRef()[patchi];
                    fvPatchScalarField& phel = hel[speciei].boundaryFieldRef()[patchi];
                    fvPatchScalarField& phvel = hvel[speciei].boundaryFieldRef()[patchi];
                    fvPatchScalarField& pTv = Tv[speciei].boundaryFieldRef()[patchi];
                    fvPatchScalarField& pzetav = zetav[speciei].boundaryFieldRef()[patchi];

                    phtMix[facei] += pY[facei]*composition().HEt(speciei, pp[facei], pTt[facei]);

                    if (pTv[facei] != 0.0)
                    {
                        if (downgradeSingleTv)
                        {
                            pTv[facei] = pTvMix[facei];
                        }

                        phv[facei] = composition().HEv(speciei, pp[facei], pTv[facei]);
                        phel[facei] = composition().HEel(speciei, pp[facei], pTv[facei]);
                        phvel[facei] = composition().HEvel(speciei, pp[facei], pTv[facei]);
                        pzetav[facei] = composition().zetav(speciei, pp[facei], pTv[facei]);

                        phvMix[facei] += pY[facei]*phv[facei];
                        phelMix[facei] += pY[facei]*phel[facei];
                        phvelMix[facei] += pY[facei]*phvel[facei];
                        phMix[facei] += pY[facei]*phvel[facei];

                        if (not downgradeSingleTv)
                        {
                            pTvMix[facei] += pX[facei]*pzetav[facei]*pTv[facei];
                        }

                        pzetavMix[facei] += pX[facei]*pzetav[facei];
                        ptotZetavMix[facei] += pX[facei]*pzetav[facei];
                    }
                }//end species loop

                phMix[facei] += phtMix[facei];

                if (not downgradeSingleTv and ptotZetavMix[facei] != 0.0)
                {
                    pTvMix[facei] /= ptotZetavMix[facei];
                }
            }//end faces loop
        }
        else
        // condition on the energy fields ... the temperatures are calculated at patches
        {
            forAll(pTt, facei)
            {
                phtMix[facei] = 0.0;
                phvMix[facei] = 0.0;
                phelMix[facei] = 0.0;
                phvelMix[facei] = 0.0;
                phMix[facei] = 0.0;
                if (not downgradeSingleTv)
                {
                    pTvMix[facei] = 0.0;
                }
                pzetavMix[facei] = 0.0;
                ptotZetavMix[facei] = 0.0;

                forAll(Y, speciei)
                {
                    const fvPatchScalarField& pY = Y[speciei].boundaryField()[patchi];
                    const fvPatchScalarField& pX = X[speciei].boundaryField()[patchi];

                    fvPatchScalarField& phv = hv[speciei].boundaryFieldRef()[patchi];
                    fvPatchScalarField& phel = hel[speciei].boundaryFieldRef()[patchi];
                    fvPatchScalarField& phvel = hvel[speciei].boundaryFieldRef()[patchi];
                    fvPatchScalarField& pTv = Tv[speciei].boundaryFieldRef()[patchi];
                    fvPatchScalarField& pzetav = zetav[speciei].boundaryFieldRef()[patchi];

                    phtMix[facei] += pY[facei]*composition().HEt(speciei, pp[facei], pTt[facei]);

                    if (pTv[facei] != 0.0)
                    {
                        if (downgradeSingleTemperature)
                        {
                            phv[facei] = composition().HEv(speciei, pp[facei], pTt[facei]);
                            phel[facei] = composition().HEel(speciei, pp[facei], pTt[facei]);
                            phvel[facei] = composition().HEvel(speciei, pp[facei], pTt[facei]);
                            pTv[facei] = pTt[facei];
                        }
                        else if (downgradeSingleTv)
                        {
                            phv[facei] = composition().HEv(speciei, pp[facei], pTvMix[facei]);
                            phel[facei] = composition().HEel(speciei, pp[facei], pTvMix[facei]);
                            phvel[facei] = composition().HEvel(speciei, pp[facei], pTvMix[facei]);
                            pTv[facei] = pTvMix[facei];
                        }
                        else if (composition().vibTempAssociativity(speciei) == -1)
                        {
                            phv[facei] = composition().HEv(speciei, pp[facei], pTv[facei]);
                            phel[facei] = composition().HEel(speciei, pp[facei], pTv[facei]);
                            phvel[facei] = composition().HEvel(speciei, pp[facei], pTv[facei]);
                            pTv[facei] = TvelEvels(speciei, phvel[facei], pp[facei], pTv[facei]);
                        }
                        else
                        {
                            const fvPatchScalarField& pTvMol = Tv[composition().vibTempAssociativity(speciei)].boundaryField()[patchi];
                            phvel[facei] = composition().HEvel(speciei, pp[facei], pTvMol[facei]);
                            pTv[facei] = pTvMol[facei];
                        }

                        phvMix[facei] += pY[facei]*phv[facei];
                        phelMix[facei] += pY[facei]*phel[facei];
                        phvelMix[facei] += pY[facei]*phvel[facei];
                        phMix[facei] += pY[facei]*phvel[facei];
                        pzetav[facei] = composition().zetav(speciei, pp[facei], pTv[facei]);

                        if (not downgradeSingleTv)
                        {
                            pTvMix[facei] += pX[facei]*pzetav[facei]*pTv[facei];
                        }

                        pzetavMix[facei] += pX[facei]*pzetav[facei];
                        ptotZetavMix[facei] += pX[facei]*pzetav[facei];
                    }
                }//end species loop

                phMix[facei] += phtMix[facei];

                if (not downgradeSingleTv and ptotZetavMix[facei] != 0.0)
                {
                    pTvMix[facei] /= ptotZetavMix[facei];
                }
            }//end faces loop
        }

        fvPatchScalarField& ppsiMix = this->psi_.boundaryFieldRef()[patchi];
        fvPatchScalarField& prhoMix = this->rho_.boundaryFieldRef()[patchi];

        //- Initialisation
        forAll(pTt, facei)
        {
            ppsiMix[facei] = 0.0;
            prhoMix[facei] = 0.0;
        }

        //- Calculation
        forAll(Y, speciei)
        {
            const fvPatchScalarField& pX = X[speciei].boundaryField()[patchi];
            const fvPatchScalarField& pTv = Tv[speciei].boundaryField()[patchi];

            forAll(pTt, facei)
            {
                if (composition().particleType(speciei) > 0)
                {
                    ppsiMix[facei] += pX[facei]*composition().psi(speciei, pp[facei], pTt[facei]);
                    prhoMix[facei] += pX[facei]*composition().rho(speciei, pp[facei], pTt[facei]);
                }
                else
                {
                    ppsiMix[facei] += pX[facei]*composition().psi(speciei, pp[facei], pTv[facei]);
                    prhoMix[facei] += pX[facei]*composition().rho(speciei, pp[facei], pTv[facei]);
                }
            }
        }
    }

    correctOverallTemperature();
}


// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -//


void Foam::rho2ReactionThermo::calculateFromDSMC
(
    const PtrList<volScalarField>& fixedY,
    const volScalarField& fixedTt,
    const PtrList<volScalarField>& fixedTv,
    const volScalarField& fixedRho,
    const word& regionName
)
{
    label zoneID = this->Tt_.mesh().cellZones().findZoneID(regionName);

    if (zoneID == -1)
    {
        FatalErrorIn("fixedFields")
            << "Cannot find region: " << regionName << nl << "in: "
            << this->Tt_.mesh().time().constant()/"polyMesh/cellZones"
            << exit(FatalError);
    }

    const cellZone& zone = this->Tt_.mesh().cellZones()[zoneID];

    if (zone.size())
    {
        PtrList<scalarList> pZone(this->Tt_.boundaryField().size());

        forAll(this->Tt_.boundaryField(), patchi)
        {
            scalarList cellsWithBoundaries(0);

            if (this->Tt_.boundaryField()[patchi].size())
            {
                const labelList& cellsNextToPatchi
                    = this->Tt_.mesh().boundary()[patchi].faceCells();

                forAll(zone, idx)
                {
                    forAll(cellsNextToPatchi, cellpatchi)
                    {
                        if (zone[idx] == cellsNextToPatchi[cellpatchi])
                        {
                            cellsWithBoundaries.append(cellpatchi);
                        }
                    }
                }
            }

            pZone.set(patchi, new scalarList(cellsWithBoundaries));
        }


        PtrList<Foam::volScalarField>& Y = composition().Y();

        scalarField& TtCells = this->Tt_.primitiveFieldRef();
        const scalarField& fixedTtCells = fixedTt.primitiveField();
        PtrList<Foam::volScalarField>& Tv = composition().Tv();
        scalarField& rhoCells = this->rho_.primitiveFieldRef();
        const scalarField& fixedRhoCells = fixedRho.primitiveField();

        scalarField& pCells = this->p_.primitiveFieldRef();
        scalarField& psiCells = this->psi_.primitiveFieldRef();
        scalarField& htCells = this->het_.primitiveFieldRef();
        PtrList<Foam::volScalarField>& hvel = composition().hevel();
        scalarField& hvelCellsMix = this->hevel().primitiveFieldRef();
        scalarField& hCells = composition().e().primitiveFieldRef();

        //- Cells values
        forAll(zone, idx)
        {
            const label celli = zone[idx];
            htCells[celli] = 0.0;
            hvelCellsMix[celli] = 0.0;
            psiCells[celli] = 0.0;

            TtCells[celli] = fixedTtCells[celli];
            rhoCells[celli] = fixedRhoCells[celli];
        }

        forAll(Y, speciei)
        {
            scalarField& YCells = Y[speciei].primitiveFieldRef();
            const scalarField& fixedYCells = fixedY[speciei].primitiveField();
            scalarField& TvCells = Tv[speciei].primitiveFieldRef();
            const scalarField& fixedTvCells = fixedTv[speciei].primitiveField();
            scalarField& hvelCells = hvel[speciei].primitiveFieldRef();

            forAll(zone, idx)
            {
                const label celli = zone[idx];
                YCells[celli] = fixedYCells[celli];
                TvCells[celli] = fixedTvCells[celli];

                htCells[celli] += YCells[celli]*composition().HEt(speciei, pCells[celli], TtCells[celli]);

                if (TvCells[celli] != 0.0)
                {
                    hvelCells[celli] = composition().HEvel(speciei, pCells[celli], TvCells[celli]);
                    hvelCellsMix[celli] += YCells[celli]*hvelCells[celli];
                }
                else
                {
                    hvelCells[celli] = 0.0;
                }

                if (composition().particleType(speciei) > 0)
                {
                    psiCells[celli] += composition().molarFraction(speciei, YCells[celli], celli)
                        *composition().psi(speciei,pCells[celli], TtCells[celli]);
                }
                else
                {
                    psiCells[celli] += composition().molarFraction(speciei, YCells[celli], celli)
                        *composition().psi(speciei,pCells[celli], TvCells[celli]);
                }
            }
        }

        forAll(zone, idx)
        {
            const label celli = zone[idx];
            hCells[celli] = htCells[celli] + hvelCellsMix[celli];
            pCells[celli] = rhoCells[celli]/psiCells[celli];
        }


        //--- boundaryField
        forAll(this->Tt_.boundaryField(), patchi)
        {
            if (this->Tt_.boundaryField()[patchi].size())
            {
                fvPatchScalarField& pTt = this->Tt_.boundaryFieldRef()[patchi];
                const fvPatchScalarField& fixedpTt = fixedTt.boundaryField()[patchi];
                fvPatchScalarField& prho = this->rho_.boundaryFieldRef()[patchi];
                const fvPatchScalarField& fixedpRho = fixedRho.boundaryField()[patchi];

                fvPatchScalarField& pp = this->p_.boundaryFieldRef()[patchi];
                fvPatchScalarField& ppsi = this->psi_.boundaryFieldRef()[patchi];
                fvPatchScalarField& pht = this->het_.boundaryFieldRef()[patchi];

                fvPatchScalarField& phvelMix = this->hevel().boundaryFieldRef()[patchi];
                fvPatchScalarField& ph = composition().e().boundaryFieldRef()[patchi];

                //- Patch values, initialisation
                forAll(pZone[patchi], idx)
                {
                    const label facei = pZone[patchi][idx];
                    pht[facei] = 0.0;
                    phvelMix[facei] = 0.0;
                    ppsi[facei] = 0.0;

                    pTt[facei] = fixedpTt[facei];
                    prho[facei] = fixedpRho[facei];
                }

                //- Patch values, calculations
                forAll(Y, speciei)
                {
                    const fvPatchScalarField& fixedpY
                        = fixedY[speciei].boundaryField()[patchi];
                    const fvPatchScalarField& fixedpTv
                        = fixedTv[speciei].boundaryField()[patchi];

                    fvPatchScalarField& pY = Y[speciei].boundaryFieldRef()[patchi];
                    fvPatchScalarField& pTv = Tv[speciei].boundaryFieldRef()[patchi];
                    fvPatchScalarField& phvel = hvel[speciei].boundaryFieldRef()[patchi];

                    forAll(pZone[patchi], idx)
                    {
                        const label facei = pZone[patchi][idx];
                        pY[facei] = fixedpY[facei];
                        pTv[facei] = fixedpTv[facei];

                        pht[facei] += pY[facei] * composition().HEt(speciei, pp[facei], pTt[facei]);

                        if (pTv[facei] != 0.0)
                        {
                            phvel[facei] = composition().HEvel(speciei,
                                pp[facei], pTv[facei]);
                            phvelMix[facei] += pY[facei] * phvel[facei];
                        }
                        else
                        {
                            phvel[facei] = 0.0;
                        }

                        if (composition().particleType(speciei) > 0)
                        {
                            ppsi[facei] += composition().molarFraction(speciei, pY[facei], facei)
                                *composition().psi(speciei, pp[facei], pTt[facei]);
                        }
                        else
                        {
                            ppsi[facei] += composition().molarFraction(speciei, pY[facei], facei)
                                *composition().psi(speciei, pp[facei], pTv[facei]);
                        }
                    }
                }//end species loop

                forAll(pZone[patchi], idx)
                {
                    const label facei = pZone[patchi][idx];
                    ph[facei] = pht[facei] + phvelMix[facei];
                    pp[facei] = prho[facei]/ppsi[facei];
                }
            }
        }//end patches loop
    }

    correctChemFractions();
}


// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -//


void Foam::rho2ReactionThermo::calculate()
{
    //- Declarations
    const bool downgradeSingleTemperature = this->downgradeSingleTemperature();
    const bool downgradeSingleVibMode = this->downgradeSingleVibMode(); // NEW VINCENT 14/03/2016 TODO ONGOING WORK
    const bool downgradeSingleTv = this->downgradeSingleTv();
    const PtrList<Foam::volScalarField>& Y = composition().Y();
    const PtrList<Foam::volScalarField>& X = composition().X();
    const scalarField& pCells = this->p_.internalField();
    const scalarField& hCellsMix = composition().e().internalField();
    const scalarField& htCellsMix = this->het_.internalField();

    PtrList<Foam::volScalarField>& Tv = composition().Tv();
    PtrList<Foam::volScalarField>& hv = composition().hev();
    PtrList<Foam::volScalarField>& hel = composition().heel();
    PtrList<Foam::volScalarField>& hvel = composition().hevel();
    PtrList<Foam::volScalarField>& zetav = composition().zetav();

    /*PtrList<PtrList<Foam::volScalarField> >& Tv_mode = composition().Tv_mode(); // NEW VINCENT 14/03/2016 TODO ONGOING WORK
    PtrList<PtrList<Foam::volScalarField> >& hevel_mode = composition().hevel_mode(); // NEW VINCENT 23/03/2016 TODO ONGOING WORK
    PtrList<PtrList<Foam::volScalarField> >& zetav_mode = composition().zetav_mode(); // NEW VINCENT 14/03/2016 TODO ONGOING WORK*/

    scalarField& TtCells = this->Tt_.primitiveFieldRef();
    scalarField& TvCellsMix = this->Tv_.primitiveFieldRef();
    scalarField& hvCellsMix = this->hevMix_.primitiveFieldRef();
    scalarField& helCellsMix = this->heelMix_.primitiveFieldRef();
    scalarField& hvelCellsMix = this->hevel().primitiveFieldRef();
    scalarField& psiCellsMix = this->psi_.primitiveFieldRef();
    scalarField& zetavCellsMix = this->zetav_.primitiveFieldRef();
    scalarField totZetavCellsMix = zetavCellsMix;


    //- Initialisations
    hvCellsMix = 0.0;
    helCellsMix = 0.0;
    if (not downgradeSingleTv)
    {
        TvCellsMix = 0.0;
        hvelCellsMix = 0.0;
    }
    zetavCellsMix = 0.0;
    totZetavCellsMix = 0.0;
    psiCellsMix = 0.0;


    //- Overall temperature calculation (single-temperature model)
    //  or trans-rotational temperature calculation (two-temperature model)
    PtrList<scalar> YList (Y.size());

    forAll(pCells, celli)
    {
        forAll(Y, speciei)
        {
            YList.set(speciei,
            new scalar(Y[speciei].internalField()[celli]));
        }

        if (downgradeSingleTemperature)
        {
            TtCells[celli] = TEs
            (
                hCellsMix[celli],
                pCells[celli],
                TtCells[celli],
                YList
            );
        }
        else if (downgradeSingleTv)
        {
            TtCells[celli] = TtEts
            (
                htCellsMix[celli],
                pCells[celli],
                TtCells[celli],
                YList
            );

            TvCellsMix[celli] = TvelEvels
            (
                hvelCellsMix[celli],
                pCells[celli],
                TvCellsMix[celli],
                YList
            );
        }
        else
        {
            TtCells[celli] = TtEts
            (
                htCellsMix[celli],
                pCells[celli],
                TtCells[celli],
                YList
            );
        }
    }//end cell loop

    forAll(this->Tt_.boundaryField(), patchi)
    {
        const volScalarField::Boundary wallPatches = this->Tt_.boundaryField();
        const fvPatchScalarField& pp = this->p_.boundaryField()[patchi];
        const fvPatchScalarField& phMix = composition().e().boundaryFieldRef()[patchi];
        fvPatchScalarField& phtMix = this->het_.boundaryFieldRef()[patchi];
        fvPatchScalarField& phvelMix = this->hevel().boundaryFieldRef()[patchi];
        fvPatchScalarField& pTt = this->Tt_.boundaryFieldRef()[patchi];
        fvPatchScalarField& pTvMix = this->Tv_.boundaryFieldRef()[patchi];

        if (pTt.fixesValue())
        // the temperature is fixed ... the energy is calculated at patches
        {
            //- Initialisation
            forAll(pTt, facei)
            {
                phtMix[facei] = 0.0;
            }

            //- Calculation
            forAll(Y, speciei)
            {
                const fvPatchScalarField& pY = Y[speciei].boundaryField()[patchi];

                forAll(pTt, facei)
                {
                    phtMix[facei] += pY[facei]*composition().HEt(speciei, pp[facei], pTt[facei]);
                }
            }
        }
        else
        // condition on the energy fields ... the temperature is calculated at patches
        {
            forAll(pTt, facei)
            {
                forAll(YList, speciei)
                {
                    YList[speciei] = Y[speciei].boundaryField()[patchi][facei];
                }

                if (downgradeSingleTemperature)
                {
                    pTt[facei] = TEs
                    (
                        phMix[facei],
                        pp[facei],
                        pTt[facei],
                        YList
                    );
                }
                else if (downgradeSingleTv)
                {
                    pTt[facei] = TtEts
                    (
                        phtMix[facei],
                        pp[facei],
                        pTt[facei],
                        YList
                    );

                    pTvMix[facei] = TvelEvels
                    (
                        phvelMix[facei],
                        pp[facei],
                        pTvMix[facei],
                        YList
                    );
                }
                else
                {
                    pTt[facei] = TtEts
                    (
                        phtMix[facei],
                        pp[facei],
                        pTt[facei],
                        YList
                    );
                }
            }//end faces loop
        }

        if (isA<wallFvPatch>(wallPatches[patchi].patch())) // NEW VINCENT 22/08/2016
        // to increase the stability of the 2T model in near-wall high Kn-number regions
        {
            forAll(pTt, facei)
            {
                if (pTt[facei] > ThighPatches_) pTt[facei] = ThighPatches_;
                else if (pTt[facei] < TlowPatches_) pTt[facei] = TlowPatches_;
            }
        }
    }//end patches loop


    //- Cells values calculation
    forAll(Y, speciei)
    {
        const scalarField& YCells = Y[speciei].internalField();
        const scalarField& XCells = X[speciei].internalField();

        scalarField& TvCells = Tv[speciei].primitiveFieldRef();
        scalarField& hvCells = hv[speciei].primitiveFieldRef();
        scalarField& helCells = hel[speciei].primitiveFieldRef();
        scalarField& hvelCells = hvel[speciei].primitiveFieldRef();
        scalarField& zetavCells = zetav[speciei].primitiveFieldRef();

        forAll(pCells, celli)
        {
            if (TvCells[celli] != 0.0)
            {
                if (downgradeSingleTemperature)
                {
                    TvCells[celli] = TtCells[celli];
                    hvelCells[celli] = composition().HEvel(speciei, pCells[celli], TvCells[celli]);
                    hvelCellsMix[celli] += YCells[celli]*hvelCells[celli];
                    zetavCells[celli] = composition().zetav(speciei, pCells[celli], TvCells[celli]);
                }
                else if (downgradeSingleTv)
                {
                    TvCells[celli] = TvCellsMix[celli];
                    hvelCells[celli] = composition().HEvel(speciei, pCells[celli], TvCells[celli]);
                    zetavCells[celli] = composition().zetav(speciei, pCells[celli], TvCells[celli]);
                }
                else if (composition().vibTempAssociativity(speciei) == -1)
                {
                    if (composition().noVibrationalTemp(speciei) == 1 or downgradeSingleVibMode)
                    {
                        if (TtCells[celli] > vibrationalCutOffTemp)
                        {
                            if (YCells[celli] > miniYforSolvingEvEqn)
                            // NEW VINCENT 08/02/2017 -> means that evEqn IS solved because Y_m large enough
                            {
                                if (sign(hvelCells[celli]) == -1) // TODO analyse and revise if necessary
                                {
                                    //Info << "hvelCells[" << celli << "]" << tab << hvelCells[celli] << endl;
                                    //Info << "B1: TvCells[" << celli << "]" << tab << TvCells[celli] << endl;

                                    hvelCells[celli] = composition().HEvel(speciei, pCells[celli], vibrationalCutOffTemp);
                                    TvCells[celli] = TvelEvels(speciei, hvelCells[celli], pCells[celli], TvCells[celli]);
                                    //Info << "A1: TvCells[" << celli << "]" << tab << TvCells[celli] << endl;
                                }
                                else
                                {
                                    //Info << "B2: TvCells[" << celli << "]" << tab << TvCells[celli] << endl;
                                    TvCells[celli] = TvelEvels(speciei, hvelCells[celli], pCells[celli], TvCells[celli]);
                                    //Info << "A2: TvCells[" << celli << "]" << tab << TvCells[celli] << endl;
                                }
                            }
                            else
                            {
                                TvCells[celli] = Tlow_; //TtCells[celli]; // TODO improve the treatment
                                //Info << "A3: TvCells[" << celli << "]" << tab << TvCells[celli] << endl;
                            }
                        }
                        else
                        {
                            TvCells[celli] = TtCells[celli];
                            //Info << "A4: TvCells[" << celli << "]" << tab << TvCells[celli] << endl;
                        }

                        hvelCellsMix[celli] += YCells[celli]*hvelCells[celli];
                        zetavCells[celli] = composition().zetav(speciei, pCells[celli], TvCells[celli]);
                    }
                    else
                    {
                        //TODO ONGOING WORK
                        /*zetavCells[celli] = 0.0;
                        TvCells[celli] = 0.0;
                        hvelCells[celli] = 0.0;

                        forAll(Tv_mode[speciei], mode)
                        {
                            scalarField& hvelmodeCells = hevel_mode[speciei][mode].internalField();
                            scalarField& TvmodeCells = Tv_mode[speciei][mode].internalField();
                            scalarField& zetavmodeCells = zetav_mode[speciei][mode].internalField();
                            TvmodeCells[celli] = composition().TvelHEvel_mode
                                (
                                    speciei, mode, hvelmodeCells[celli],
                                    pCells[celli], TvmodeCells[celli]
                                );
                            zetavmodeCells[celli] = composition().zetav_mode(speciei, mode, pCells[celli], TvmodeCells[celli]);
                            zetavCells[celli] += zetavmodeCells[celli];
                            hvelCells[celli] += hvelmodeCells[celli];
                            TvCells[celli] += TvmodeCells[celli]*zetavmodeCells[celli];
                        }

                        TvCells[celli] /= zetavCells[celli];*/
                    }
                }
                else
                {
                    // If no molecule in the flow-field, then necessarily the 1-T solver is run using
                    // downgradeSingleTemperature because the electronic energy equation is not
                    // implemented and it cannot be substituted by a vibrational energy equation.
                    // Therefore, the following always works. VINCENT 08/08/2016
                    TvCells[celli] = Tv[composition().vibTempAssociativity(speciei)].internalField()[celli];
                    hvelCells[celli] = composition().HEvel(speciei, pCells[celli], TvCells[celli]);
                    zetavCells[celli] = composition().zetav(speciei, pCells[celli], TvCells[celli]);
                    hvelCellsMix[celli] += YCells[celli]*hvelCells[celli];
                }

                hvCells[celli] = composition().HEv(speciei, pCells[celli], TvCells[celli]);
                helCells[celli] = composition().HEel(speciei, pCells[celli], TvCells[celli]);

                hvCellsMix[celli] += YCells[celli]*hvCells[celli];
                helCellsMix[celli] += YCells[celli]*helCells[celli];
                if (not downgradeSingleTv)
                {
                    TvCellsMix[celli] += XCells[celli]*zetavCells[celli]*TvCells[celli];
                }
                zetavCellsMix[celli] += XCells[celli]*zetavCells[celli];
                totZetavCellsMix[celli] += XCells[celli]*zetavCells[celli];
            }

            if (composition().particleType(speciei) > 0)
            {
                psiCellsMix[celli] += XCells[celli]*composition().psi(speciei, pCells[celli], TtCells[celli]);
            }
            else
            {
                psiCellsMix[celli] += XCells[celli]*composition().psi(speciei, pCells[celli], TvCells[celli]);
            }
        }//end cells loop
    }//end species loop

    forAll(pCells, celli)
    {
        if (not downgradeSingleTv and totZetavCellsMix[celli] != 0.0)
        {
            TvCellsMix[celli] /= totZetavCellsMix[celli];
        }
    }

    //- Patch values calculation
    forAll(this->Tt_.boundaryField(), patchi)
    {
        const fvPatchScalarField& pp = this->p_.boundaryField()[patchi];
        const fvPatchScalarField& pTt = this->Tt_.boundaryField()[patchi];

        fvPatchScalarField& phMix = composition().e().boundaryFieldRef()[patchi];
        fvPatchScalarField& phtMix = this->het_.boundaryFieldRef()[patchi];
        fvPatchScalarField& phvMix = this->hevMix_.boundaryFieldRef()[patchi];
        fvPatchScalarField& phelMix = this->heelMix_.boundaryFieldRef()[patchi];
        fvPatchScalarField& phvelMix = this->hevel().boundaryFieldRef()[patchi];
        fvPatchScalarField& pTvMix = this->Tv_.boundaryFieldRef()[patchi];
        fvPatchScalarField& pzetavMix = this->zetav_.boundaryFieldRef()[patchi];
        fvPatchScalarField ptotZetavMix = pzetavMix;

        fvPatchScalarField& ppsiMix = this->psi_.boundaryFieldRef()[patchi];

        //- Re-set of mixture quantities boundaryField
        forAll(pTt, facei)
        {
            phvMix[facei] = 0.0;
            phelMix[facei] = 0.0;
            if (not downgradeSingleTv)
            {
                phvelMix[facei] = 0.0;
                pTvMix[facei] = 0.0;
            }
            pzetavMix[facei] = 0.0;
            ptotZetavMix[facei] = 0.0;

            ppsiMix[facei] = 0.0;
        }

        if (pTt.fixesValue())
        // the temperature is fixed ... the energy is calculated at patches
        {
            forAll(Y, speciei)
            {
                const fvPatchScalarField& pY = Y[speciei].boundaryField()[patchi];
                const fvPatchScalarField& pX = X[speciei].boundaryField()[patchi];
                fvPatchScalarField& pTv = Tv[speciei].boundaryFieldRef()[patchi];
                fvPatchScalarField& phv = hv[speciei].boundaryFieldRef()[patchi];
                fvPatchScalarField& phel = hel[speciei].boundaryFieldRef()[patchi];
                fvPatchScalarField& phvel = hvel[speciei].boundaryFieldRef()[patchi];
                fvPatchScalarField& pzetav = zetav[speciei].boundaryFieldRef()[patchi];

                forAll(pTt, facei)
                {
                    if (pTv[facei] != 0.0)
                    {
                        if (downgradeSingleTemperature)
                        {
                            phv[facei] = composition().HEv(speciei, pp[facei], pTt[facei]);
                            phel[facei] = composition().HEel(speciei, pp[facei], pTt[facei]);
                            phvel[facei] = composition().HEvel(speciei, pp[facei], pTt[facei]);
                            phvelMix[facei] += pY[facei]*phvel[facei];
                            pzetav[facei] = composition().zetav(speciei, pp[facei], pTv[facei]);
                            pTvMix[facei] += pX[facei]*pzetav[facei]*pTv[facei];
                        }
                        else if (downgradeSingleTv)
                        {
                            pTv[facei] = pTvMix[facei];
                            phv[facei] = composition().HEv(speciei, pp[facei], pTv[facei]);
                            phel[facei] = composition().HEel(speciei, pp[facei], pTv[facei]);
                            phvel[facei] = composition().HEvel(speciei, pp[facei], pTv[facei]);
                            pzetav[facei] = composition().zetav(speciei, pp[facei], pTv[facei]);
                        }
                        else
                        {
                            phvel[facei] = composition().HEvel(speciei, pp[facei], pTv[facei]);
                            phv[facei] = composition().HEv(speciei, pp[facei], pTv[facei]);
                            phel[facei] = composition().HEel(speciei, pp[facei], pTv[facei]);
                            phvelMix[facei] += pY[facei]*phvel[facei];
                            pzetav[facei] = composition().zetav(speciei, pp[facei], pTv[facei]);
                            pTvMix[facei] += pX[facei]*pzetav[facei]*pTv[facei];
                        }

                        phvMix[facei] += pY[facei]*phv[facei];
                        phelMix[facei] += pY[facei]*phel[facei];
                        pzetavMix[facei] += pX[facei]*pzetav[facei];
                        ptotZetavMix[facei] += pX[facei]*pzetav[facei];
                    }
                }//end faces loop
            }//end species loop

            forAll(pTt, facei)
            {
                phMix[facei] = phtMix[facei] + phvelMix[facei];

                if (not downgradeSingleTv and ptotZetavMix[facei] != 0.0)
                {
                    pTvMix[facei] /= ptotZetavMix[facei];
                }
            }
        }
        else
        // condition on the energy fields ... the temperatures are calculated at patches
        {
            forAll(pTt, facei)
            {
                forAll(Y, speciei)
                {
                    const fvPatchScalarField& pY = Y[speciei].boundaryField()[patchi];
                    const fvPatchScalarField& pX = X[speciei].boundaryField()[patchi];
                    const fvPatchScalarField& phv = hv[speciei].boundaryField()[patchi];
                    const fvPatchScalarField& phel = hel[speciei].boundaryField()[patchi];

                    fvPatchScalarField& phvel = hvel[speciei].boundaryFieldRef()[patchi];
                    fvPatchScalarField& pTv = Tv[speciei].boundaryFieldRef()[patchi];
                    fvPatchScalarField& pzetav = zetav[speciei].boundaryFieldRef()[patchi];

                    if (pTv[facei] != 0.0)
                    {
                        if (downgradeSingleTemperature)
                        {
                            pTv[facei] = pTt[facei];
                            phvel[facei] = composition().HEvel(speciei, pp[facei], pTv[facei]);
                            pzetav[facei] = composition().zetav(speciei, pp[facei], pTv[facei]);
                            phvelMix[facei] += pY[facei]*phvel[facei];
                            pTvMix[facei] += pX[facei]*pzetav[facei]*pTv[facei];
                        }
                        else if (downgradeSingleTv)
                        {
                            pTv[facei] = pTvMix[facei];
                            phvel[facei] = composition().HEvel(speciei, pp[facei], pTv[facei]);
                            pzetav[facei] = composition().zetav(speciei, pp[facei], pTv[facei]);
                        }
                        else if (composition().vibTempAssociativity(speciei) == -1)
                        {
                            if (composition().noVibrationalTemp(speciei) == 1 or downgradeSingleVibMode)
                            {
                                if (pTt[facei] > vibrationalCutOffTemp)
                                {
                                    if (pY[facei] > miniYforSolvingEvEqn)
                                    // NEW VINCENT 08/02/2017 -> < means that evEqn IS solved because Y_m large enough
                                    {
                                        if (sign(phvel[facei]) == -1) // NEW VINCENT 09/02/2017
                                        {
                                            phvel[facei] = composition().HEvel(speciei, pp[facei], vibrationalCutOffTemp);
                                        }

                                        pTv[facei] = TvelEvels(speciei, phvel[facei], pp[facei], pTv[facei]);
                                    }
                                    else
                                    {
                                        pTv[facei] = TlowPatches_; //pTt[facei]; // TODO improve the treatment
                                    }
                                }

                                pzetav[facei] = composition().zetav(speciei, pp[facei], pTv[facei]);
                                phvelMix[facei] += pY[facei]*phvel[facei];
                                pTvMix[facei] += pX[facei]*pzetav[facei]*pTv[facei];
                            }
                            else
                            {
                                //TODO ONGOING WORK
                                /*pzetav[facei] = 0.0;
                                pTv[facei] = 0.0;

                                forAll(Tv_mode[speciei], mode)
                                {
                                    Tv_mode[speciei][mode].boundaryField()[patchi][facei] =
                                        composition().TvelHEvel_mode(speciei, mode, phvel[facei], pp[facei], Tv_mode[speciei][mode].boundaryField()[patchi][facei]);
                                    zetav_mode[speciei][mode].boundaryField()[patchi][facei] =
                                        composition().zetav_mode(speciei, mode, pp[facei], Tv_mode[speciei][mode].boundaryField()[patchi][facei]);
                                    pzetav[facei] += zetav_mode[speciei][mode].boundaryField()[patchi][facei];
                                    pTv[facei] += Tv_mode[speciei][mode].internalField()[facei]*zetav_mode[speciei][mode].boundaryField()[patchi][facei];
                                }

                                pTv[facei] /= pzetav[facei];*/
                            }
                        }
                        else
                        // ions, electrons, and atoms if Eel is on
                        {
                            const fvPatchScalarField& pTvMol = Tv[composition().vibTempAssociativity(speciei)].boundaryField()[patchi];
                            pTv[facei] = pTvMol[facei];
                            pzetav[facei] = composition().zetav(speciei, pp[facei], pTv[facei]);
                        }

                        phvMix[facei] += pY[facei]*phv[facei];
                        phelMix[facei] += pY[facei]*phel[facei];
                        pzetavMix[facei] += pX[facei]*pzetav[facei];
                        ptotZetavMix[facei] += pX[facei]*pzetav[facei];
                    }
                }//end species loop

                if (not downgradeSingleTv and ptotZetavMix[facei] != 0.0)
                {
                    pTvMix[facei] /= ptotZetavMix[facei];
                }
            }//end faces loop
        }

        forAll(Y, speciei)
        {
            const fvPatchScalarField& pX = X[speciei].boundaryField()[patchi];
            const fvPatchScalarField& pTv = Tv[speciei].boundaryField()[patchi];

            forAll(pTt, facei)
            {
                if (composition().particleType(speciei) > 0)
                {
                    ppsiMix[facei] += pX[facei]*composition().psi(speciei, pp[facei], pTt[facei]);
                }
                else
                {
                    ppsiMix[facei] += pX[facei]*composition().psi(speciei, pp[facei], pTv[facei]);
                }
            }//end faces loop
        }//end species loop
    }//end patches loop
}


// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -//


void Foam::rho2ReactionThermo::calculateLight()
{
    //- Declarations
    const bool downgradeSingleTemperature = this->downgradeSingleTemperature_;
    const bool downgradeSingleTv = this->downgradeSingleTv_;
    const PtrList<Foam::volScalarField>& Y = composition().Y();
    const PtrList<Foam::volScalarField>& X = composition().X();
    const scalarField& pCells = this->p_.internalField();
    const scalarField& hCellsMix = composition().e().internalField();
    const scalarField& htCellsMix = this->het_.internalField();

    PtrList<Foam::volScalarField>& Tv = composition().Tv();
    PtrList<Foam::volScalarField>& hvel = composition().hevel();

    scalarField& TtCells = this->Tt_.primitiveFieldRef();
    scalarField& TvCellsMix = this->Tv_.primitiveFieldRef();
    scalarField& hvelCellsMix = this->hevel().primitiveFieldRef();
    scalarField& psiCellsMix = this->psi_.primitiveFieldRef();

    // Overall temperature calculation (single-temperature model)
    // or trans-rotational temperature calculation (two-temperature model)
    PtrList<scalar> YList(Y.size());

    forAll(pCells, celli)
    {
        forAll(Y, speciei)
        {
            YList.set
            (
                speciei,
                new scalar(Y[speciei].internalField()[celli])
            );
        }

        if (downgradeSingleTemperature)
        {
            TtCells[celli] = TEs
            (
                hCellsMix[celli],
                pCells[celli],
                TtCells[celli],
                YList
            );
        }
        else if (downgradeSingleTv)
        {
            TtCells[celli] = TtEts
            (
                htCellsMix[celli],
                pCells[celli],
                TtCells[celli],
                YList
            );

            TvCellsMix[celli] = TvelEvels
            (
                hvelCellsMix[celli],
                pCells[celli],
                TvCellsMix[celli],
                YList
            );
        }
        else
        {
            TtCells[celli] = TtEts
            (
                htCellsMix[celli],
                pCells[celli],
                TtCells[celli],
                YList
            );
        }
    }//end cells loop

    forAll(this->Tt_.boundaryField(), patchi)
    {
        const volScalarField::Boundary wallPatches = this->Tt_.boundaryField(); // NEW VINCENT 23/08/2016
        const fvPatchScalarField& pp = this->p_.boundaryField()[patchi];
        const fvPatchScalarField& phMix = composition().e().boundaryField()[patchi];
        const fvPatchScalarField& phvelMix = this->hevel().boundaryField()[patchi]; // NEW VINCENT 16/08/2016
        fvPatchScalarField& phtMix = this->het_.boundaryFieldRef()[patchi];
        fvPatchScalarField& pTt = this->Tt_.boundaryFieldRef()[patchi];
        fvPatchScalarField& pTvMix = this->Tv_.boundaryFieldRef()[patchi]; // NEW VINCENT 16/08/2016

        if (pTt.fixesValue())
        // the temperature is fixed ... the energy is calculated at patches
        {
            //- Initialisation
            forAll(pTt, facei)
            {
                phtMix[facei] = 0.0;
            }

            //- Calculation
            forAll(Y, speciei)
            {
                const fvPatchScalarField& pY = Y[speciei].boundaryField()[patchi];

                forAll(pTt, facei)
                {
                    phtMix[facei] += pY[facei]*composition().HEt(speciei, pp[facei], pTt[facei]);
                }
            }
        }
        else
        // condition on the energy fields ... the temperature is calculated at patches
        {
            forAll(pTt, facei)
            {
                forAll(Y, speciei)
                {
                    YList[speciei] = Y[speciei].boundaryField()[patchi][facei];
                }

                if (downgradeSingleTemperature)
                {
                    pTt[facei] = TEs
                    (
                        phMix[facei],
                        pp[facei],
                        pTt[facei],
                        YList
                    );
                }
                else if (downgradeSingleTv) // NEW VINCENT 16/08/2016
                {
                    pTt[facei] = TtEts
                    (
                        phtMix[facei],
                        pp[facei],
                        pTt[facei],
                        YList
                    );

                    pTvMix[facei] = TvelEvels
                    (
                        phvelMix[facei],
                        pp[facei],
                        pTvMix[facei],
                        YList
                    );
                }
                else
                {
                    pTt[facei] = TtEts
                    (
                        phtMix[facei],
                        pp[facei],
                        pTt[facei],
                        YList
                    );
                }
            }//end faces loop
        }

        if (isA<wallFvPatch>(wallPatches[patchi].patch())) // NEW VINCENT 22/08/2016
        // to increase the stability of the 2T model in near-wall high-Kn number regions
        {
            forAll(pTt, facei)
            {
                if (pTt[facei] > ThighPatches_) pTt[facei] = ThighPatches_;
                else if (pTt[facei] < TlowPatches_) pTt[facei] = TlowPatches_;
            }
        }
    }//end patches loop


    //- Re-set of mixture quantities internalField
    if (not downgradeSingleTv) hvelCellsMix = 0.0;
    psiCellsMix = 0.0;

    //- Cells values calculation
    forAll(Y, speciei)
    {
        const scalarField& YCells = Y[speciei].internalField();
        const scalarField& XCells = X[speciei].internalField();
        scalarField& TvCells = Tv[speciei].primitiveFieldRef();
        scalarField& hvelCells = hvel[speciei].primitiveFieldRef();

        forAll(YCells, celli)
        {
            if (TvCells[celli] != 0.0)
            {
                if (downgradeSingleTemperature)
                {
                    TvCells[celli] = TtCells[celli];
                    hvelCells[celli] = composition().HEvel(speciei, pCells[celli], TvCells[celli]);
                    hvelCellsMix[celli] += YCells[celli]*hvelCells[celli];
                }
                else if (downgradeSingleTv)
                {
                    TvCells[celli] = TvCellsMix[celli];
                    hvelCells[celli] = composition().HEvel(speciei, pCells[celli], TvCells[celli]);
                }
                else if (composition().vibTempAssociativity(speciei) == -1) // NEW VINCENT 11/08/2016
                {
                    if (TtCells[celli] > vibrationalCutOffTemp)
                    {
                        if (YCells[celli] > miniYforSolvingEvEqn)
                        // NEW VINCENT 08/02/2017 -> < means that evEqn IS solved because Y_m large enough
                        {
                            if (sign(hvelCells[celli]) == -1)
                            {
                                hvelCells[celli] = composition().HEvel(speciei, pCells[celli], vibrationalCutOffTemp); // NEW VINCENT 09/02/2017
                            }

                            TvCells[celli] = TvelEvels(speciei, hvelCells[celli], pCells[celli], TvCells[celli]);
                        }
                        else
                        {
                            TvCells[celli] = Tlow_; //TtCells[celli]; // TODO improve the treatment
                        }
                    }
                    else
                    {
                        TvCells[celli] = TtCells[celli];
                    }

                    hvelCellsMix[celli] += YCells[celli]*hvelCells[celli];
                }
                else
                {
                    // If no molecule in the flow-field, then necessarily the 1-T solver is run using
                    // downgradeSingleTemperature because the electronic energy equation is not
                    // implemented and it cannot be substituted by a vibrational energy equation.
                    // Therefore, the following always works. VINCENT 08/08/2016
                    TvCells[celli] = Tv[composition().vibTempAssociativity(speciei)].internalField()[celli];
                    hvelCells[celli] = composition().HEvel(speciei, pCells[celli], TvCells[celli]);
                    hvelCellsMix[celli] += YCells[celli]*hvelCells[celli];
                }
            }

            if (composition().particleType(speciei) > 0)
            {
                psiCellsMix[celli] += XCells[celli]*composition().psi(speciei, pCells[celli], TtCells[celli]);
            }
            else
            {
                psiCellsMix[celli] += XCells[celli]*composition().psi(speciei, pCells[celli], TvCells[celli]);
            }
        }//end cells loop
    }//end species loop


    //- Patch values calculation
    forAll(this->Tt_.boundaryField(), patchi)
    {
        const fvPatchScalarField& pp = this->p_.boundaryField()[patchi];
        const fvPatchScalarField& pTt = this->Tt_.boundaryField()[patchi];
        const fvPatchScalarField& pTvMix = this->Tv_.boundaryField()[patchi]; // NEW VINCENT 15/08/2016
        fvPatchScalarField& phtMix = this->het_.boundaryFieldRef()[patchi];
        fvPatchScalarField& phMix = composition().e().boundaryFieldRef()[patchi];
        fvPatchScalarField& phvelMix = this->hevel().boundaryFieldRef()[patchi];
        fvPatchScalarField& ppsiMix = this->psi_.boundaryFieldRef()[patchi];

        //- Re-set of mixture quantities boundaryField
        forAll(pTt, facei)
        {
            if (not downgradeSingleTv) phvelMix[facei] = 0.0;
            ppsiMix[facei] = 0.0;
        }

        if (pTt.fixesValue())
        // the temperature is fixed ... the energy is calculated at patches
        {
            forAll(Y, speciei)
            {
                const fvPatchScalarField& pY = Y[speciei].boundaryField()[patchi];
                fvPatchScalarField& pTv = Tv[speciei].boundaryFieldRef()[patchi];
                fvPatchScalarField& phvel = hvel[speciei].boundaryFieldRef()[patchi];

                forAll(pTt, facei)
                {
                    if (pTv[facei] != 0.0)
                    {
                        if (downgradeSingleTemperature)
                        {
                            phvel[facei] = composition().HEvel(speciei, pp[facei], pTt[facei]);
                            phvelMix[facei] += pY[facei]*phvel[facei];
                        }
                        else if (downgradeSingleTv)
                        {
                            pTv[facei] = pTvMix[facei];
                            phvel[facei] = composition().HEvel(speciei, pp[facei], pTv[facei]);
                        }
                        else
                        {
                            phvel[facei] = composition().HEvel(speciei, pp[facei], pTv[facei]);
                            phvelMix[facei] += pY[facei]*phvel[facei];
                        }
                    }
                }//end faces loop
            }//end species loop

            forAll(pTt, facei)
            {
                phMix[facei] = phtMix[facei] + phvelMix[facei];
            }
        }
        else
        // condition on the energy field ... the temperature is calculated at patches
        {
            forAll(Y, speciei)
            {
                const fvPatchScalarField& pY = Y[speciei].boundaryField()[patchi];
                fvPatchScalarField& phvel = hvel[speciei].boundaryFieldRef()[patchi];
                fvPatchScalarField& pTv = Tv[speciei].boundaryFieldRef()[patchi];

                forAll(pTt, facei)
                {
                    if (pTv[facei] != 0.0)
                    {
                        if (downgradeSingleTemperature)
                        {
                            pTv[facei] = pTt[facei];
                            phvel[facei] = composition().HEvel(speciei, pp[facei], pTv[facei]);
                            phvelMix[facei] += pY[facei]*phvel[facei];
                        }
                        else if (downgradeSingleTv)
                        {
                            pTv[facei] = pTvMix[facei];
                            phvel[facei] = composition().HEvel(speciei, pp[facei], pTv[facei]);
                        }
                        else if (composition().vibTempAssociativity(speciei) == -1) // NEW VINCENT 11/08/2016
                        {
                            if (pTt[facei] > vibrationalCutOffTemp)
                            {
                                if (pY[facei] > miniYforSolvingEvEqn)
                                // NEW VINCENT 08/02/2017 -> < means that evEqn IS solved because Y_m large enough
                                {
                                    if (sign(phvel[facei]) == -1) // NEW VINCENT 09/02/2017
                                    {
                                        phvel[facei] = composition().HEvel(speciei, pp[facei], vibrationalCutOffTemp);
                                    }

                                    pTv[facei] = TvelEvels(speciei, phvel[facei], pp[facei], pTv[facei]);
                                }
                                else
                                {
                                    pTv[facei] = TlowPatches_; //pTt[facei]; // TODO improve the treatment
                                }
                            }
                            else
                            {
                                pTv[facei] = pTt[facei];
                            }

                            phvelMix[facei] += pY[facei]*phvel[facei];
                        }
                        else
                        {
                            pTv[facei] = Tv[composition().vibTempAssociativity(speciei)].boundaryField()[patchi][facei];
                            phvel[facei] = composition().HEvel(speciei, pp[facei], pTv[facei]);
                            phvelMix[facei] += pY[facei]*phvel[facei];
                        }
                    }
                }//end faces loop
            }//end species loop

            forAll(pTt, facei)
            {
                phtMix[facei] = phMix[facei] - phvelMix[facei]; // NEW VINCENT 25/09/2016
            }
        }

        forAll(Y, speciei)
        {
            const fvPatchScalarField& pX = X[speciei].boundaryField()[patchi];
            const fvPatchScalarField& pTv = Tv[speciei].boundaryField()[patchi];

            forAll(pTt, facei)
            {
                if (composition().particleType(speciei) > 0)
                {
                    ppsiMix[facei] += pX[facei]*composition().psi(speciei, pp[facei], pTt[facei]);
                }
                else
                {
                    ppsiMix[facei] += pX[facei]*composition().psi(speciei, pp[facei], pTv[facei]);
                }
            }//end faces loop
        }//end species loop
    }//end patches loop
}


// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -//


void Foam::rho2ReactionThermo::calculateHEVK()
{
    // Tve is set according to the vibTempAssociativity table for particles whose vib energy is not solved

    //- Declaration
    const bool downgradeSingleTv = this->downgradeSingleTv();

    const scalarField& pCells = this->p_.internalField();
    const scalarField& TvCellsMix = this->Tv_.internalField();
    PtrList<Foam::volScalarField>& Tv = composition().Tv();
    PtrList<Foam::volScalarField>& hvel = composition().hevel();

    //- Cells values
    forAll(Tv, speciei)
    {
        scalarField& TvCells = Tv[speciei].primitiveFieldRef();
        scalarField& hvelCells = hvel[speciei].primitiveFieldRef();

        if (TvCells[0] != 0.0)
        {
            if (downgradeSingleTv)
            {
                forAll(pCells, celli)
                {
                    TvCells[celli] = TvCellsMix[celli];
                    hvelCells[celli] = composition().HEvel(speciei, pCells[celli], TvCells[celli]);
                }
            }
            else if (composition().vibTempAssociativity(speciei) != -1) // NEW VINCENT 04/08/2016
            {
                forAll(pCells, celli)
                {
                    TvCells[celli] = Tv[composition().vibTempAssociativity(speciei)].internalField()[celli];
                    hvelCells[celli] = composition().HEvel(speciei, pCells[celli], TvCells[celli]);
                }
            }
        }

        //- Patch values
        forAll(this->p_.boundaryField(), patchi)
        {
            fvPatchScalarField& pTv = Tv[speciei].boundaryFieldRef()[patchi];
            const fvPatchScalarField& pTvMix = this->Tv_.boundaryField()[patchi];

            if (pTv.size() != 0) // NEW VINCENT 04/08/2016
            {
                if (pTv[0] != 0.0)
                {
                    const fvPatchScalarField& pp = this->p_.boundaryField()[patchi];
                    fvPatchScalarField& phvel = hvel[speciei].boundaryFieldRef()[patchi];

                    if (downgradeSingleTv)
                    {
                        forAll(pTv, facei)
                        {
                            pTv[facei] = pTvMix[facei];
                            phvel[facei] = composition().HEvel(speciei, pp[facei], pTv[facei]);
                        }
                    }
                    else if (composition().vibTempAssociativity(speciei) != -1) // NEW VINCENT 04/08/2016
                    {
                        forAll(pTv, facei)
                        {
                            pTv[facei] = Tv[composition().vibTempAssociativity(speciei)].boundaryField()[patchi][facei];
                            phvel[facei] = composition().HEvel(speciei, pp[facei], pTv[facei]);
                        }
                    }
                }
            }
        }//end patches loop
    }//end species loop
}


void Foam::rho2ReactionThermo::limitVelocityAtWallBoundary
(
    volVectorField::Boundary& Ubdry
)
{
    //- Increase the stability of the 2T model in near-wall high Kn-number
    //  regions
    
    const volScalarField::Boundary wallPatches = this->Tt_.boundaryField();

    //- Patch values
    forAll(Ubdry, patchi)
    {
        if (isA<wallFvPatch>(wallPatches[patchi].patch()))
        {
            fvPatchVectorField& pU = Ubdry[patchi];

            forAll(pU, facei)
            {
                if (pU[facei].component(0) > UhighPatches_)
                {
                    pU[facei].component(0) = UhighPatches_;
                }
                if (pU[facei].component(1) > UhighPatches_)
                {
                    pU[facei].component(1) = UhighPatches_;
                }
                if (pU[facei].component(2) > UhighPatches_)
                {
                    pU[facei].component(2) = UhighPatches_;
                }
            }
        }
    }//end patches loop
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::rho2ReactionThermo::rho2ReactionThermo
(
    const fvMesh& mesh,
    const word& phaseName
)
:
    rho2Thermo(mesh, phaseName),
    partialThermoName_
    (
        transportToTypedef(word(subDict("thermoType").lookup("transport")))
    )
{
    if (this->downgradeSingleTv_)
    {
        hevelMix_ = new volScalarField
        (
            IOobject
            (
                this->phasePropertyName("hevel"),
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimEnergy/dimMass,
            this->hevelMix2BoundaryTypes(),
            this->hevelMix2BoundaryBaseTypes()
        );
    }

    if (isDict("temperatureBounds"))
    {
        Tlow_ =
            subDict("temperatureBounds").lookupOrDefault<scalar>("Tlow", 100.0);
        Thigh_ =
            subDict("temperatureBounds").lookupOrDefault<scalar>
            (
                "Thigh",
                40000.0
            );
        TlowPatches_ =
            subDict("temperatureBounds").lookupOrDefault<scalar>
            (
                "TlowPatches",
                Tlow_
            );
        ThighPatches_ =
            subDict("temperatureBounds").lookupOrDefault<scalar>
            (
                "ThighPatches",
                Thigh_
            );
    }
    else
    {
        Tlow_ = 100;
        Thigh_ = 40000;
        TlowPatches_ = Tlow_;
        ThighPatches_ = Thigh_;
    }

    if (isDict("velocityBounds"))
    {
        UhighPatches_ =
            subDict("velocityBounds").lookupOrDefault<scalar>
            (
                "UhighPatches",
                Foam::GREAT
            );
    }
    else
    {
        UhighPatches_ = Foam::GREAT;
    }
}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::rho2ReactionThermo> Foam::rho2ReactionThermo::New
(
    const fvMesh& mesh,
    const word& phaseName
)
{
    return basic2Thermo::New<rho2ReactionThermo>(mesh, phaseName);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::rho2ReactionThermo::~rho2ReactionThermo()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::rho2ReactionThermo& Foam::rho2ReactionThermo::lookup2ReactionThermo
(
    const fvPatchScalarField& pf
)
{
    return pf.db().lookupObject<rho2ReactionThermo>(dictName);
}


void Foam::rho2ReactionThermo::correctFractions()
{
    if (debug)
    {
        Info<< "entering rho2ReactionThermo::correctFractions()" << endl;
    }

    correctChemFractions();

    if (debug)
    {
        Info<< "exiting rho2ReactionThermo::correctFractions()" << endl;
    }
}


void Foam::rho2ReactionThermo::correctFromDSMC
(
    const PtrList<volScalarField>& fixedY,
    const volScalarField& fixedTt,
    const PtrList<volScalarField>& fixedTv,
    const volScalarField& fixedRho,
    const word& regionName
)
{
    calculateFromDSMC
    (
        fixedY,
        fixedTt,
        fixedTv,
        fixedRho,
        regionName
    );
}


void Foam::rho2ReactionThermo::correct2T()
{
    if (debug)
    {
        Info<< "entering rho2ReactionThermo::correct2T()" << endl;
    }

    calculate();

    if (debug)
    {
        Info<< "exiting rho2ReactionThermo::correct2T()" << endl;
    }
}


void Foam::rho2ReactionThermo::correct2T_Light()
{
    if (debug)
    {
        Info<< "entering rho2ReactionThermo::correct2T_Light()" << endl;
    }

    calculateLight();

    if (debug)
    {
        Info<< "exiting rho2ReactionThermo::correct2T_Light()" << endl;
    }
}


void Foam::rho2ReactionThermo::correctHEVK()
{
    if (debug)
    {
        Info<< "entering rho2ReactionThermo::correctHEVK()" << endl;
    }

    calculateHEVK();

    if (debug)
    {
        Info<< "exiting rho2ReactionThermo::correctHEVK()" << endl;
    }
}


void Foam::rho2ReactionThermo::initialise2T()
{
    if (debug)
    {
        Info<< "entering rho2ReactionThermo::initialise2T()" << endl;
    }

    initialise();

    if (debug)
    {
        Info<< "exiting rho2ReactionThermo::initialise2T()" << endl;
    }
}


void Foam::rho2ReactionThermo::initialise2T_Light()
{
    if (debug)
    {
        Info<< "entering rho2ReactionThermo::initialise2T_Light()" << endl;
    }

    initialiseLight();

    if (debug)
    {
        Info<< "exiting rho2ReactionThermo::initialise2T_Light()" << endl;
    }
}


Foam::word Foam::rho2ReactionThermo::partialThermoName()
{
    return partialThermoName_;
}


// ************************************************************************* //
