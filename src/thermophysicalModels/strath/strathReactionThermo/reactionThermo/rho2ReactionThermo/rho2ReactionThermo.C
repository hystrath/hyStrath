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


void Foam::rho2ReactionThermo::correctVolChemFractions()
{
    const scalar Runi = constant::physicoChemical::R.value();
    const scalarField& rhoCells = this->rho_.internalField();
    
    PtrList<Foam::volScalarField>& Y = composition().Y();
    PtrList<Foam::volScalarField>& X = composition().X();
    PtrList<Foam::volScalarField>& nD = composition().nD();
    PtrList<Foam::volScalarField>& pD = composition().pD();

    volScalarField& Wmix = composition().Wmix();
    volScalarField& Rmix = this->RMix_;
    
    scalarField& WmixCells = Wmix.primitiveFieldRef();
    scalarField& RmixCells = Rmix.primitiveFieldRef();
    
    scalarField YtotCells(rhoCells.size(), 0.0);

    //- Bound Y between 0 and 1 and compute Wmix from solved Y
    WmixCells = 0.0;
    
    forAll(Y, speciei)
    {
        const scalar W = composition().W(speciei);
        scalarField& YCells = Y[speciei].primitiveFieldRef();
        
        forAll(YCells, celli)
        {
            YCells[celli] = min(max(YCells[celli], 0.0), 1.0);
            YtotCells[celli] += YCells[celli];
            WmixCells[celli] += YCells[celli]/W;
        }
    }
    
    WmixCells = 1.0/WmixCells;
    
    if (composition().contains("e-"))
    {
        const scalarField& XeCells = composition().X("e-");
        
        scalarField XtotCells(rhoCells.size(), 0.0);
        scalarField XtotIonCells(rhoCells.size(), 0.0);
        
        //- Compute X from known Y and Wmix
        forAll(Y, speciei)
        {
            const scalar W = composition().W(speciei);
            const scalarField& YCells = Y[speciei].internalField();
            
            scalarField& XCells = X[speciei].primitiveFieldRef();
            
            forAll(YCells, celli)
            {
                XCells[celli] = YCells[celli]*WmixCells[celli]/W;
                
                if (composition().isNeutral(speciei))
                {
                    XtotCells[celli] += XCells[celli];
                }
                else if (composition().isIon(speciei))
                {
                    XtotIonCells[celli] += XCells[celli];
                }
            }
        }
        
        //- Compute X_corrected ensuring that
        //    the sum of all X equals 1
        //    charge conservation
        WmixCells = 0.0;
        
        forAll(Y, speciei)
        {
            const scalar W = composition().W(speciei);
            scalarField& XCells = X[speciei].primitiveFieldRef();
            
            forAll(XCells, celli)
            {
                if (composition().isNeutral(speciei))
                {
                    XCells[celli] *= (1.0-2.0*XeCells[celli])/XtotCells[celli];
                }
                else if
                (
                    composition().isIon(speciei) and XtotIonCells[celli] > SMALL
                )
                {
                    XCells[celli] *= XeCells[celli]/XtotIonCells[celli];
                }
                
                //- Compute Wmix_corrected from X_corrected
                WmixCells[celli] += XCells[celli]*W;
            }
        }
        
        RmixCells = 1000.0*Runi/WmixCells;
        
        forAll(Y, speciei)
        {
            const scalar W = composition().W(speciei);
            const scalarField& XCells = X[speciei].internalField();
            scalarField& YCells = Y[speciei].primitiveFieldRef();
            scalarField& nDCells = nD[speciei].primitiveFieldRef();
            scalarField& pDCells = pD[speciei].primitiveFieldRef();
            
            forAll(YCells, celli)
            {
                //- Compute Y_corrected from X_corrected and Wmix_corrected
                YCells[celli] = XCells[celli]*W/WmixCells[celli];
                
                //- Compute nD and pD from Y_corrected and rho
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
            }
        }
    }
    else
    {
        WmixCells = 0.0;
        
        //- Compute Y_corrected ensuring that the sum of all Y equals 1
        forAll(Y, speciei)
        {
            const scalar W = composition().W(speciei);
            scalarField& YCells = Y[speciei].primitiveFieldRef();
            scalarField& nDCells = nD[speciei].primitiveFieldRef();
            scalarField& pDCells = pD[speciei].primitiveFieldRef();
            
            forAll(YCells, celli)
            {
                YCells[celli] /= YtotCells[celli];
                
                //- Compute Wmix_corrected from Y_corrected
                WmixCells[celli] += YCells[celli]/W;
                
                //- Compute nD and pD from Y_corrected and rho
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
            }
        }
        
        RmixCells = 1000.0*Runi*WmixCells;
        WmixCells = 1.0/WmixCells;
        
        //- Compute X from Y_corrected
        forAll(Y, speciei)
        {
            const scalar W = composition().W(speciei);
            const scalarField& YCells = Y[speciei].internalField();
            scalarField& XCells = X[speciei].primitiveFieldRef();
            
            forAll(YCells, celli)
            {
                XCells[celli] = YCells[celli]*WmixCells[celli]/W;
            }
        }
    }
}


void Foam::rho2ReactionThermo::correctBdrChemFractions()
{
    const PtrList<Foam::volScalarField>& Y = composition().Y();
    PtrList<Foam::volScalarField>& X = composition().X();
    PtrList<Foam::volScalarField>& nD = composition().nD();
    PtrList<Foam::volScalarField>& pD = composition().pD();
    volScalarField& Wmix = composition().Wmix();
    volScalarField& Rmix = this->RMix_;
    
    //- Boundary patches
    forAll(this->rho_.boundaryField(), patchi)
    {
        const fvPatchScalarField& prho = this->rho_.boundaryField()[patchi];
        
        fvPatchScalarField& pWmix = Wmix.boundaryFieldRef()[patchi];
        fvPatchScalarField& pRmix = Rmix.boundaryFieldRef()[patchi];

        //- Initialisation
        pWmix = 0.0;

        forAll(Y, speciei)
        {
            const scalar W = composition().W(speciei);
            const fvPatchScalarField& pY = Y[speciei].boundaryField()[patchi];
            fvPatchScalarField& pnD = nD[speciei].boundaryFieldRef()[patchi];
            fvPatchScalarField& ppD = pD[speciei].boundaryFieldRef()[patchi];

            forAll(pY, facei)
            {
                pWmix[facei] += pY[facei]/W;
                
                pnD[facei] =
                    composition().numberDensity
                    (
                        speciei,
                        pY[facei],
                        prho[facei]
                    );
                    
                ppD[facei] =
                    composition().partialDensity(pY[facei], prho[facei]);
            }//end faces loop
        }//end species loop
        
        pRmix = 1000.0*constant::physicoChemical::R.value()*pWmix;
        pWmix = 1.0/pWmix;
        
        forAll(Y, speciei)
        {
            const fvPatchScalarField& pY = Y[speciei].boundaryField()[patchi];
            fvPatchScalarField& pX = X[speciei].boundaryFieldRef()[patchi];
            
            //- Compute X from Y and Wmix
            forAll(pY, facei)
            {
                pX[facei] =
                    composition().molarFraction
                    (
                        speciei,
                        pY[facei],
                        patchi,
                        facei
                    );
            }//end faces loop
        }//end species loop
        
    }//end patches loop
}


void Foam::rho2ReactionThermo::correctOverallTemperature()
{
    //- Declaration
    const PtrList<Foam::volScalarField>& Y = composition().Y();
    const PtrList<Foam::volScalarField>& X = composition().X();
    const PtrList<Foam::volScalarField>& Tv = composition().Tv();
    const scalarField& pCells = this->p().internalField();
    const scalarField& TtCells = this->T().internalField();

    PtrList<Foam::volScalarField>& zetar = composition().zetar();
    PtrList<Foam::volScalarField>& zetav = composition().zetav();
    PtrList<Foam::volScalarField>& zetael = composition().zetael();

    volScalarField& Tov = this->Tov();
    volScalarField& zetarMix = this->zetar();
    volScalarField& zetaelMix = this->zetael();
    
    volScalarField numTovMix = this->Tov();
    volScalarField denTovMix = this->Tov();

    scalarField& zetarCellsMix = zetarMix.primitiveFieldRef();
    scalarField& zetaelCellsMix = zetaelMix.primitiveFieldRef();
    scalarField& numTovCellsMix = numTovMix.primitiveFieldRef();
    scalarField& denTovCellsMix = denTovMix.primitiveFieldRef();
    scalarField& TovCells = Tov.primitiveFieldRef();

    //- Initialisation
    zetarCellsMix = 0.0;
    zetaelCellsMix = 0.0;
    numTovCellsMix = 0.0;
    denTovCellsMix = 0.0;

    forAll(this->p().boundaryField(), patchi)
    {
        fvPatchScalarField& pzetarMix =
            this->zetar().boundaryFieldRef()[patchi];
        fvPatchScalarField& pzetaelMix =
            this->zetael().boundaryFieldRef()[patchi];
        fvPatchScalarField& pnumTovMix = numTovMix.boundaryFieldRef()[patchi];
        fvPatchScalarField& pdenTovMix = denTovMix.boundaryFieldRef()[patchi];

        pzetarMix = 0.0;
        pzetaelMix = 0.0;
        pnumTovMix = 0.0;
        pdenTovMix = 0.0;
    }

    //- Update degrees of freedom for each energy mode
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

            if (composition().isElectron(speciei))
            {
                numTovCellsMix[celli] += 3.0*TvCells[celli]*YCells[celli];
                denTovCellsMix[celli] += 3.0*YCells[celli];
            }
            else
            {
                numTovCellsMix[celli] +=
                    (
                        (3.0 + zetarCells[celli])*TtCells[celli]
                      + (zetavCells[celli] + zetaelCells[celli])*TvCells[celli]
                    )*YCells[celli];
                denTovCellsMix[celli] +=
                    (
                        3.0 + zetarCells[celli] + zetavCells[celli]
                      + zetaelCells[celli]
                    )*YCells[celli];
            }
        }//end cells loop

        //- Patch values calculation
        forAll(this->p().boundaryField(), patchi)
        {
            const fvPatchScalarField& pp = this->p().boundaryField()[patchi];
            const fvPatchScalarField& pTt = this->T().boundaryField()[patchi];
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
            fvPatchScalarField& pnumTovMix = numTovMix.boundaryFieldRef()[patchi];
            fvPatchScalarField& pdenTovMix = denTovMix.boundaryFieldRef()[patchi];

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
                    pnumTovMix[facei] += 3.0*pTv[facei]*pY[facei];
                    pdenTovMix[facei] += 3.0*pY[facei];
                }
                else
                {
                    pnumTovMix[facei] +=
                        (
                            (3.0 + pzetar[facei])*pTt[facei] 
                          + (pzetav[facei] + pzetael[facei])*pTv[facei]
                        )*pY[facei];
                    pdenTovMix[facei] +=
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
        TovCells[celli] = numTovCellsMix[celli]/denTovCellsMix[celli];
    }

    forAll(this->p().boundaryField(), patchi)
    {
        const fvPatchScalarField& pnumTovMix = numTovMix.boundaryField()[patchi];
        const fvPatchScalarField& pdenTovMix = denTovMix.boundaryField()[patchi];
        fvPatchScalarField& pTov = this->Tov_.boundaryFieldRef()[patchi];

        forAll(pTov, facei)
        {
            pTov[facei] = pnumTovMix[facei]/pdenTovMix[facei];
        }
    }
}


// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //


void Foam::rho2ReactionThermo::initialise()
{
    //- Declarations
    const bool downgradeSingleTemperature = this->downgradeSingleTemperature();
    const bool downgradeSingleTv = this->downgradeSingleTv();
    const PtrList<Foam::volScalarField>& Y = composition().Y();
    const PtrList<Foam::volScalarField>& X = composition().X();
    const scalarField& pCells = this->p().internalField();
    const scalarField& TtCells = this->T().internalField();
    const scalarField& RmixCells = this->RMix_.primitiveFieldRef();

    PtrList<Foam::volScalarField>& hv = composition().hev();
    PtrList<Foam::volScalarField>& hel = composition().heel();
    PtrList<Foam::volScalarField>& hvel = composition().hevel();
    PtrList<Foam::volScalarField>& Tv = composition().Tv();
    PtrList<Foam::volScalarField>& zetav = composition().zetav();
    PtrList<Foam::volScalarField>& pP = composition().pP();
    PtrList<Foam::volScalarField>& Cvtr = composition().Cvtr();
    PtrList<Foam::volScalarField>& Cvvel = composition().Cvvel();

    // ABORTIVE WORK
    //PtrList<PtrList<Foam::volScalarField> >& hvel_mode = composition().hevel_mode();
    //PtrList<PtrList<Foam::volScalarField> >& Tv_mode = composition().Tv_mode();
    //PtrList<PtrList<Foam::volScalarField> >& zetav_mode = composition().zetav_mode();

    scalarField& TvCellsMix = this->Tv().primitiveFieldRef();
    scalarField& htCellsMix = this->het().primitiveFieldRef();
    scalarField& hvCellsMix = this->hevMix_.primitiveFieldRef();
    scalarField& helCellsMix = this->heelMix_.primitiveFieldRef();
    scalarField& hvelCellsMix = this->hevel().primitiveFieldRef();
    scalarField& hCellsMix = composition().e().primitiveFieldRef();
    scalarField& zetavCellsMix = this->zetav_.primitiveFieldRef();
    scalarField totZetavCellsMix = zetavCellsMix;
    scalarField& psiCellsMix = this->psi_.primitiveFieldRef();
    scalarField& rhoCellsMix = this->rho_.primitiveFieldRef();
    scalarField& CvtrCellsMix = this->CvtrMix_.primitiveFieldRef();
    scalarField& CvvelCellsMix = this->CvvelMix_.primitiveFieldRef();
    scalarField& CvCellsMix = this->CvMix_.primitiveFieldRef();
    scalarField& CptrCellsMix = this->CptrMix_.primitiveFieldRef();
    scalarField& CpvelCellsMix = this->CpvelMix_.primitiveFieldRef();
    scalarField& CpCellsMix = this->CpMix_.primitiveFieldRef();

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
        CvtrCellsMix[celli] = 0.0;
        CvvelCellsMix[celli] = 0.0;
        CptrCellsMix[celli] = 0.0;
        CpvelCellsMix[celli] = 0.0;
        
        //- Calculation
        forAll(Y, speciei)
        {
            const scalarField& YCells = Y[speciei].internalField();
            htCellsMix[celli] += YCells[celli]
               *composition().HEt(speciei, pCells[celli], TtCells[celli]);
        }

        forAll(Y, speciei)
        {
            const scalar R = composition().R(speciei);
            const scalarField& YCells = Y[speciei].internalField();
            const scalarField& XCells = X[speciei].internalField();

            scalarField& hvCells = hv[speciei].primitiveFieldRef();
            scalarField& helCells = hel[speciei].primitiveFieldRef();
            scalarField& hvelCells = hvel[speciei].primitiveFieldRef();
            scalarField& TvCells = Tv[speciei].primitiveFieldRef();
            scalarField& pPCells = pP[speciei].primitiveFieldRef();
            scalarField& zetavCells = zetav[speciei].primitiveFieldRef();
            scalarField& CvtrCells = Cvtr[speciei].primitiveFieldRef();
            scalarField& CvvelCells = Cvvel[speciei].primitiveFieldRef();

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

                hvCells[celli] =
                    composition().HEv
                    (
                        speciei,
                        pCells[celli],
                        TvCells[celli]
                    );
                helCells[celli] =
                    composition().HEel
                    (
                        speciei,
                        pCells[celli],
                        TvCells[celli]
                    );
                hvelCells[celli] =
                    composition().HEvel
                    (
                        speciei,
                        pCells[celli],
                        TvCells[celli]
                    );
                // ABORTIVE WORK
                /*forAll(hvel_mode[speciei], vibMode)
                {
                    hvel_mode[speciei][vibMode].internalField()[celli] =
                        composition().HEvel_mode
                        (
                            speciei,
                            vibMode,
                            pCells[celli],
                            Tv_mode[speciei][vibMode].internalField()[celli]
                        );
                }*/

                zetavCells[celli] =
                    composition().zetav
                    (
                        speciei,
                        pCells[celli],
                        TvCells[celli]
                    );

                hvCellsMix[celli] += YCells[celli]*hvCells[celli];
                helCellsMix[celli] += YCells[celli]*helCells[celli];
                hvelCellsMix[celli] += YCells[celli]*hvelCells[celli];
                hCellsMix[celli] += YCells[celli]*hvelCells[celli];
                if (not downgradeSingleTv)
                {
                    TvCellsMix[celli] +=
                        XCells[celli]*zetavCells[celli]*TvCells[celli];
                }
                zetavCellsMix[celli] += XCells[celli]*zetavCells[celli];
                totZetavCellsMix[celli] += XCells[celli]*zetavCells[celli];
            }

            CvtrCells[celli] =
                composition().Cv_t(speciei, pCells[celli], TtCells[celli]);
            CvvelCells[celli] =
                composition().Cv_vel(speciei, pCells[celli], TvCells[celli]);
            CvtrCellsMix[celli] += YCells[celli]*CvtrCells[celli];
            CvvelCellsMix[celli] += YCells[celli]*CvvelCells[celli];
            CptrCellsMix[celli] += YCells[celli]*CvtrCells[celli];
            CpvelCellsMix[celli] += YCells[celli]*CvvelCells[celli];
                
            scalar pPToRhomix = 0.0;

            if (composition().isHeavySpecies(speciei))
            {
                CptrCellsMix[celli] += YCells[celli]*R;

                pPToRhomix = YCells[celli]
                    /composition().psi(speciei, pCells[celli], TtCells[celli]);
            }
            else
            {
                CpvelCellsMix[celli] += YCells[celli]*R;
                
                pPToRhomix = YCells[celli]
                    /composition().psi(speciei, pCells[celli], TvCells[celli]);
            }

            pPCells[celli] = pPToRhomix;
            psiCellsMix[celli] += pPToRhomix;
        }//end species loop
        
        psiCellsMix[celli] = 1.0/psiCellsMix[celli];
        rhoCellsMix[celli] = pCells[celli]*psiCellsMix[celli];
        
        CvCellsMix[celli] = CvtrCellsMix[celli] + CvvelCellsMix[celli];
        CpCellsMix[celli] = CvCellsMix[celli] + RmixCells[celli];

        if (not downgradeSingleTv and totZetavCellsMix[celli] != 0.0)
        {
            TvCellsMix[celli] /= totZetavCellsMix[celli];
        }

        hCellsMix[celli] += htCellsMix[celli];
    }//end cells loop
    
    //- Volume calculation for partial pressures
    forAll(Y, speciei)
    {
        scalarField& pPCells = pP[speciei].primitiveFieldRef();
        forAll(pCells, celli)
        {
            pPCells[celli] *= rhoCellsMix[celli];
        }
    }

    //- Patch values
    forAll(this->T().boundaryField(), patchi)
    {
        const fvPatchScalarField& pp = this->p().boundaryField()[patchi];
        const fvPatchScalarField& pTt = this->T().boundaryField()[patchi];
        const fvPatchScalarField& pRmix = this->RMix_.boundaryField()[patchi];

        fvPatchScalarField& phtMix = this->het().boundaryFieldRef()[patchi];
        fvPatchScalarField& phvMix = this->hevMix_.boundaryFieldRef()[patchi];
        fvPatchScalarField& phelMix = this->heelMix_.boundaryFieldRef()[patchi];
        fvPatchScalarField& phvelMix = this->hevel().boundaryFieldRef()[patchi];
        fvPatchScalarField& phMix =
            composition().e().boundaryFieldRef()[patchi];
        fvPatchScalarField& pTvMix = this->Tv().boundaryFieldRef()[patchi];
        fvPatchScalarField& pzetavMix = this->zetav_.boundaryFieldRef()[patchi];
        fvPatchScalarField ptotZetavMix = pzetavMix;

        if (pTt.fixesValue())
        {
            // the temperature is fixed ... the energy is calculated at patches
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
                    const fvPatchScalarField& pY =
                        Y[speciei].boundaryField()[patchi];
                    const fvPatchScalarField& pX =
                        X[speciei].boundaryField()[patchi];

                    fvPatchScalarField& phv =
                        hv[speciei].boundaryFieldRef()[patchi];
                    fvPatchScalarField& phel =
                        hel[speciei].boundaryFieldRef()[patchi];
                    fvPatchScalarField& phvel =
                        hvel[speciei].boundaryFieldRef()[patchi];
                    fvPatchScalarField& pTv =
                        Tv[speciei].boundaryFieldRef()[patchi];
                    fvPatchScalarField& pzetav =
                        zetav[speciei].boundaryFieldRef()[patchi];

                    phtMix[facei] += pY[facei]
                        *composition().HEt(speciei, pp[facei], pTt[facei]);

                    if (pTv[facei] != 0.0)
                    {
                        if (downgradeSingleTv)
                        {
                            pTv[facei] = pTvMix[facei];
                        }

                        phv[facei] =
                            composition().HEv(speciei, pp[facei], pTv[facei]);
                        phel[facei] =
                            composition().HEel(speciei, pp[facei], pTv[facei]);
                        phvel[facei] =
                            composition().HEvel(speciei, pp[facei], pTv[facei]);
                        pzetav[facei] =
                            composition().zetav(speciei, pp[facei], pTv[facei]);

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
        {
            // condition on the energy fields ... the temperatures are
            // calculated at patches
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
                    const fvPatchScalarField& pY =
                        Y[speciei].boundaryField()[patchi];
                    const fvPatchScalarField&
                        pX = X[speciei].boundaryField()[patchi];

                    fvPatchScalarField& phv =
                        hv[speciei].boundaryFieldRef()[patchi];
                    fvPatchScalarField& phel =
                        hel[speciei].boundaryFieldRef()[patchi];
                    fvPatchScalarField& phvel =
                        hvel[speciei].boundaryFieldRef()[patchi];
                    fvPatchScalarField& pTv =
                        Tv[speciei].boundaryFieldRef()[patchi];
                    fvPatchScalarField& pzetav =
                        zetav[speciei].boundaryFieldRef()[patchi];

                    phtMix[facei] += pY[facei]
                        *composition().HEt(speciei, pp[facei], pTt[facei]);

                    if (pTv[facei] > SMALL)
                    {
                        if (downgradeSingleTemperature)
                        {
                            phv[facei] =
                                composition().HEv
                                (
                                    speciei,
                                    pp[facei],
                                    pTt[facei]
                                );
                            phvel[facei] =
                                composition().HEvel
                                (
                                    speciei,
                                    pp[facei],
                                    pTt[facei]
                                );
                            pTv[facei] = pTt[facei];
                        }
                        else if (downgradeSingleTv)
                        {
                            phv[facei] =
                                composition().HEv
                                (
                                    speciei,
                                    pp[facei],
                                    pTvMix[facei]
                                );
                            phel[facei] =
                                composition().HEel
                                (
                                    speciei,
                                    pp[facei],
                                    pTvMix[facei]
                                );
                            phvel[facei] =
                                composition().HEvel
                                (
                                    speciei,
                                    pp[facei],
                                    pTvMix[facei]
                                );
                            pTv[facei] = pTvMix[facei];
                        }
                        else if (composition().vibTempAssociativity(speciei) == -1)
                        {
                            phv[facei] =
                                composition().HEv
                                (
                                    speciei,
                                    pp[facei],
                                    pTv[facei]
                                );
                            phvel[facei] =
                                composition().HEvel
                                (
                                    speciei,
                                    pp[facei],
                                    pTv[facei]
                                );
                            pTv[facei] =
                                TvelEvels
                                (
                                    speciei,
                                    phvel[facei],
                                    pp[facei],
                                    pTv[facei]
                                );
                        }
                        else
                        {
                            const fvPatchScalarField& pTvMol =
                                Tv[composition().vibTempAssociativity(speciei)]
                                    .boundaryField()[patchi];
                            phvel[facei] =
                                composition().HEvel
                                (
                                    speciei,
                                    pp[facei],
                                    pTvMol[facei]
                                );
                            pTv[facei] = pTvMol[facei];
                        }

                        phvMix[facei] += pY[facei]*phv[facei];
                        phelMix[facei] += pY[facei]*phel[facei];
                        phvelMix[facei] += pY[facei]*phvel[facei];
                        phMix[facei] += pY[facei]*phvel[facei];
                        pzetav[facei] =
                            composition().zetav(speciei, pp[facei], pTv[facei]);

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
        fvPatchScalarField& pCvtrMix =
            this->CvtrMix_.boundaryFieldRef()[patchi];
        fvPatchScalarField& pCvvelMix =
            this->CvvelMix_.boundaryFieldRef()[patchi];
        fvPatchScalarField& pCptrMix =
            this->CptrMix_.boundaryFieldRef()[patchi];
        fvPatchScalarField& pCpvelMix =
            this->CpvelMix_.boundaryFieldRef()[patchi];    
        fvPatchScalarField& pCvMix = this->CvMix_.boundaryFieldRef()[patchi];
        fvPatchScalarField& pCpMix = this->CpMix_.boundaryFieldRef()[patchi];

        //- Initialisation
        forAll(pTt, facei)
        {
            ppsiMix[facei] = 0.0;
            pCvtrMix[facei] = 0.0;
            pCvvelMix[facei] = 0.0;
            pCptrMix[facei] = 0.0;
            pCpvelMix[facei] = 0.0;
        }

        //- Calculation
        forAll(Y, speciei)
        {
            const scalar R = composition().R(speciei);
            const fvPatchScalarField& pY = Y[speciei].boundaryField()[patchi];
            const fvPatchScalarField& pTv = Tv[speciei].boundaryField()[patchi];
            
            fvPatchScalarField& ppP = pP[speciei].boundaryFieldRef()[patchi];
            fvPatchScalarField& pCvtr =
                Cvtr[speciei].boundaryFieldRef()[patchi];
            fvPatchScalarField& pCvvel =
                Cvvel[speciei].boundaryFieldRef()[patchi];

            forAll(pTt, facei)
            {
                pCvtr[facei] =
                    composition().Cv_t(speciei, pp[facei], pTt[facei]);
                pCvvel[facei] =
                    composition().Cv_vel(speciei, pp[facei], pTv[facei]);
                pCvtrMix[facei] += pY[facei]*pCvtr[facei];
                pCvvelMix[facei] += pY[facei]*pCvvel[facei];
                pCptrMix[facei] += pY[facei]*pCvtr[facei];
                pCpvelMix[facei] += pY[facei]*pCvvel[facei];
                
                scalar pPToRhoMix = 0.0;
                
                if (composition().isHeavySpecies(speciei))
                {
                    pCptrMix[facei] += pY[facei]*R;
                    
                    pPToRhoMix = pY[facei]
                        /composition().psi(speciei, pp[facei], pTt[facei]);
                }
                else
                {
                    pCpvelMix[facei] += pY[facei]*R;
                    
                    pPToRhoMix = pY[facei]
                        /composition().psi(speciei, pp[facei], pTv[facei]);
                }
                
                ppP[facei] = pPToRhoMix;
                ppsiMix[facei] += pPToRhoMix;
                
            }//end faces loop
        }//end species loop
        
        ppsiMix = 1.0/ppsiMix;
        prhoMix = pp*ppsiMix;
        
        forAll(Y, speciei)
        {
            fvPatchScalarField& ppP = pP[speciei].boundaryFieldRef()[patchi];
            forAll(ppP, facei)
            {
                ppP[facei] *= prhoMix[facei];
            }
        }
        
        pCvMix = pCvtrMix + pCvvelMix;
        pCpMix = pCvMix + pRmix;
        
    }//end patches loop
    
    // Calculation of the electron pressure 
    if (composition().contains("e-"))
    {
        this->pe_ = composition().pP("e-");
    }

    correctVolChemFractions();
    correctBdrChemFractions();
    
    correctOverallTemperature();
}


// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //


void Foam::rho2ReactionThermo::calculateFromDSMC
(
    const PtrList<volScalarField>& fixedY,
    const volScalarField& fixedTt,
    const PtrList<volScalarField>& fixedTv,
    const volScalarField& fixedRho,
    const word& regionName
)
{
    label zoneID = this->T().mesh().cellZones().findZoneID(regionName);

    if (zoneID == -1)
    {
        FatalErrorIn("fixedFields")
            << "Cannot find region: " << regionName << nl << "in: "
            << this->T().mesh().time().constant()/"polyMesh/cellZones"
            << exit(FatalError);
    }

    const cellZone& zone = this->T().mesh().cellZones()[zoneID];

    if (zone.size())
    {
        PtrList<scalarList> pZone(this->T().boundaryField().size());

        forAll(this->T().boundaryField(), patchi)
        {
            scalarList cellsWithBoundaries(0);

            if (this->T().boundaryField()[patchi].size())
            {
                const labelList& cellsNextToPatchi =
                    this->T().mesh().boundary()[patchi].faceCells();

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

        scalarField& TtCells = this->T().primitiveFieldRef();
        const scalarField& fixedTtCells = fixedTt.primitiveField();
        PtrList<Foam::volScalarField>& Tv = composition().Tv();
        scalarField& rhoCells = this->rho_.primitiveFieldRef();
        const scalarField& fixedRhoCells = fixedRho.primitiveField();

        scalarField& pCells = this->p().primitiveFieldRef();
        scalarField& psiCells = this->psi_.primitiveFieldRef();
        scalarField& htCells = this->het().primitiveFieldRef();
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

                htCells[celli] += YCells[celli]
                    *composition().HEt(speciei, pCells[celli], TtCells[celli]);

                if (TvCells[celli] > SMALL)
                {
                    hvelCells[celli] =
                        composition().HEvel
                        (
                            speciei,
                            pCells[celli],
                            TvCells[celli]
                        );
                    hvelCellsMix[celli] += YCells[celli]*hvelCells[celli];
                }
                else
                {
                    hvelCells[celli] = 0.0;
                }

                if (composition().isHeavySpecies(speciei))
                {
                    psiCells[celli] += YCells[celli]
                      / composition().psi
                        (
                            speciei,
                            pCells[celli],
                            TtCells[celli]
                        );
                }
                else
                {
                    psiCells[celli] += YCells[celli]
                      / composition().psi
                        (
                            speciei,
                            pCells[celli],
                            TvCells[celli]
                        );
                }
            }// end zone loop
        }// end species loop

        forAll(zone, idx)
        {
            const label celli = zone[idx];
            
            psiCells[celli] = 1.0/psiCells[celli];
            pCells[celli] = rhoCells[celli]/psiCells[celli];
            hCells[celli] = htCells[celli] + hvelCellsMix[celli];
        }


        //--- boundaryField
        forAll(this->T().boundaryField(), patchi)
        {
            if (this->T().boundaryField()[patchi].size())
            {
                fvPatchScalarField& pTt = this->T().boundaryFieldRef()[patchi];
                const fvPatchScalarField& fixedpTt =
                    fixedTt.boundaryField()[patchi];
                fvPatchScalarField& prho =
                    this->rho_.boundaryFieldRef()[patchi];
                const fvPatchScalarField& fixedpRho =
                    fixedRho.boundaryField()[patchi];

                fvPatchScalarField& pp = this->p().boundaryFieldRef()[patchi];
                fvPatchScalarField& ppsi =
                    this->psi_.boundaryFieldRef()[patchi];
                fvPatchScalarField& pht =
                    this->het().boundaryFieldRef()[patchi];

                fvPatchScalarField& phvelMix =
                    this->hevel().boundaryFieldRef()[patchi];
                fvPatchScalarField& ph =
                    composition().e().boundaryFieldRef()[patchi];

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

                    fvPatchScalarField& pY =
                        Y[speciei].boundaryFieldRef()[patchi];
                    fvPatchScalarField& pTv =
                        Tv[speciei].boundaryFieldRef()[patchi];
                    fvPatchScalarField& phvel =
                        hvel[speciei].boundaryFieldRef()[patchi];

                    forAll(pZone[patchi], idx)
                    {
                        const label facei = pZone[patchi][idx];
                        pY[facei] = fixedpY[facei];
                        pTv[facei] = fixedpTv[facei];

                        pht[facei] += pY[facei]
                            *composition().HEt(speciei, pp[facei], pTt[facei]);

                        if (pTv[facei] > SMALL)
                        {
                            phvel[facei] =
                                composition().HEvel
                                (
                                    speciei,
                                    pp[facei],
                                    pTv[facei]
                                );
                            phvelMix[facei] += pY[facei] * phvel[facei];
                        }
                        else
                        {
                            phvel[facei] = 0.0;
                        }

                        if (composition().isHeavySpecies(speciei))
                        {
                            ppsi[facei] += pY[facei]
                              / composition().psi(speciei, pp[facei], pTt[facei]);
                        }
                        else
                        {
                            ppsi[facei] += pY[facei]
                              / composition().psi(speciei, pp[facei], pTv[facei]);
                        }
                    }
                }//end species loop

                forAll(pZone[patchi], idx)
                {
                    const label facei = pZone[patchi][idx];
                    
                    ppsi[facei] = 1.0/ppsi[facei];
                    pp[facei] = prho[facei]/ppsi[facei];
                    ph[facei] = pht[facei] + phvelMix[facei];
                }
            }
        }//end patches loop
    }
    
    correctVolChemFractions();
    correctBdrChemFractions();
}


// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //


void Foam::rho2ReactionThermo::calculate()
{
    //- Declarations
    const bool downgradeSingleTemperature = this->downgradeSingleTemperature();
    const bool downgradeSingleVibMode = this->downgradeSingleVibMode();
    const bool downgradeSingleTv = this->downgradeSingleTv();
    const PtrList<Foam::volScalarField>& Y = composition().Y();
    const PtrList<Foam::volScalarField>& X = composition().X();
    const scalarField& pCells = this->p().internalField();
    const scalarField& rhoCells = this->rho_.internalField();
    const scalarField& hCellsMix = composition().e().internalField();
    const scalarField& htCellsMix = this->het().internalField();
    const scalarField& RmixCells = this->RMix_.primitiveFieldRef();

    PtrList<Foam::volScalarField>& Tv = composition().Tv();
    PtrList<Foam::volScalarField>& hv = composition().hev();
    PtrList<Foam::volScalarField>& hel = composition().heel();
    PtrList<Foam::volScalarField>& hvel = composition().hevel();
    PtrList<Foam::volScalarField>& zetav = composition().zetav();
    PtrList<Foam::volScalarField>& pP = composition().pP();
    PtrList<Foam::volScalarField>& Cvtr = composition().Cvtr();
    PtrList<Foam::volScalarField>& Cvvel = composition().Cvvel();

    // ABORTIVE WORK
    /*PtrList<PtrList<Foam::volScalarField> >& Tv_mode = composition().Tv_mode();
    PtrList<PtrList<Foam::volScalarField> >& hevel_mode = composition().hevel_mode();
    PtrList<PtrList<Foam::volScalarField> >& zetav_mode = composition().zetav_mode();*/

    scalarField& TtCells = this->T().primitiveFieldRef();
    scalarField& TvCellsMix = this->Tv().primitiveFieldRef();
    scalarField& hvCellsMix = this->hevMix_.primitiveFieldRef();
    scalarField& helCellsMix = this->heelMix_.primitiveFieldRef();
    scalarField& hvelCellsMix = this->hevel().primitiveFieldRef();
    scalarField& psiCellsMix = this->psi_.primitiveFieldRef();
    scalarField& zetavCellsMix = this->zetav_.primitiveFieldRef();
    scalarField totZetavCellsMix = zetavCellsMix;
    scalarField& CvtrCellsMix = this->CvtrMix_.primitiveFieldRef();
    scalarField& CvvelCellsMix = this->CvvelMix_.primitiveFieldRef();
    scalarField& CvCellsMix = this->CvMix_.primitiveFieldRef();
    scalarField& CptrCellsMix = this->CptrMix_.primitiveFieldRef();
    scalarField& CpvelCellsMix = this->CpvelMix_.primitiveFieldRef();
    scalarField& CpCellsMix = this->CpMix_.primitiveFieldRef();


    //- Initialisation
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
    CvtrCellsMix = 0.0;
    CvvelCellsMix = 0.0;
    CptrCellsMix = 0.0;
    CpvelCellsMix = 0.0;


    //- Overall temperature calculation (single-temperature model)
    //  or trans-rotational temperature calculation (two-temperature model)
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
            TtCells[celli] =
                TEs
                (
                    hCellsMix[celli],
                    pCells[celli],
                    TtCells[celli],
                    YList
                );
        }
        else
        {
            TtCells[celli] = 
                TtEts
                (
                    htCellsMix[celli],
                    pCells[celli],
                    TtCells[celli],
                    YList
                );

            if (downgradeSingleTv)
            {
                TvCellsMix[celli] =
                    TvelEvels
                    (
                        hvelCellsMix[celli],
                        pCells[celli],
                        TvCellsMix[celli],
                        YList
                    );
            }
        }
    }//end cells loop

    forAll(this->T().boundaryField(), patchi)
    {
        const volScalarField::Boundary wallPatches = this->T().boundaryField();
        const fvPatchScalarField& pp = this->p().boundaryField()[patchi];
        const fvPatchScalarField& phMix =
            composition().e().boundaryField()[patchi];
        const fvPatchScalarField& phtMix = this->het().boundaryField()[patchi];
        const fvPatchScalarField& phvelMix =
            this->hevel().boundaryField()[patchi];
        fvPatchScalarField& pTt = this->T().boundaryFieldRef()[patchi];
        fvPatchScalarField& pTvMix = this->Tv().boundaryFieldRef()[patchi];

        if (pTt.fixesValue())
        {
            // The temperature is fixed ... the energy is calculated at patches
            // using the fixedTREnergy boundary condition for et via the
            // == operator. The updates occurs when 
            // et.correctBoundaryConditions() is called at solver level
        }
        else
        {
            // condition on the energy fields ... the temperature is calculated
            // at patches
            forAll(pTt, facei)
            {
                forAll(Y, speciei)
                {
                    YList[speciei] = Y[speciei].boundaryField()[patchi][facei];
                }

                if (downgradeSingleTemperature)
                {
                    pTt[facei] =
                        TEs
                        (
                            phMix[facei],
                            pp[facei],
                            pTt[facei],
                            YList
                        );
                }
                else
                {
                    pTt[facei] =
                        TtEts
                        (
                            phtMix[facei],
                            pp[facei],
                            pTt[facei],
                            YList
                        );

                    if (downgradeSingleTv)
                    {
                        pTvMix[facei] =
                            TvelEvels
                            (
                                phvelMix[facei],
                                pp[facei],
                                pTvMix[facei],
                                YList
                            );
                    }
                }
            }//end faces loop
        }

        if (isA<wallFvPatch>(wallPatches[patchi].patch()))
        {
            // to increase the stability of the 2T model in near-wall high
            // Kn-number regions
            forAll(pTt, facei)
            {
                if (pTt[facei] > ThighPatches_)
                {
                    pTt[facei] = ThighPatches_;
                }
                else if (pTt[facei] < TlowPatches_)
                {
                    pTt[facei] = TlowPatches_;
                }
            }
        }
    }//end patches loop


    //- Cells values calculation
    forAll(Y, speciei)
    {
        const scalar R = composition().R(speciei);
        const scalarField& YCells = Y[speciei].internalField();
        const scalarField& XCells = X[speciei].internalField();

        scalarField& TvCells = Tv[speciei].primitiveFieldRef();
        scalarField& hvCells = hv[speciei].primitiveFieldRef();
        scalarField& helCells = hel[speciei].primitiveFieldRef();
        scalarField& hvelCells = hvel[speciei].primitiveFieldRef();
        scalarField& zetavCells = zetav[speciei].primitiveFieldRef();
        scalarField& pPCells = pP[speciei].primitiveFieldRef();
        scalarField& CvtrCells = Cvtr[speciei].primitiveFieldRef();
        scalarField& CvvelCells = Cvvel[speciei].primitiveFieldRef();

        forAll(pCells, celli)
        {
            if (TvCells[celli] > SMALL)
            {
                if (downgradeSingleTemperature)
                {
                    TvCells[celli] = TtCells[celli];
                    hvelCells[celli] =
                        composition().HEvel
                        (
                            speciei,
                            pCells[celli],
                            TvCells[celli]
                        );
                    hvelCellsMix[celli] += YCells[celli]*hvelCells[celli];
                    zetavCells[celli] =
                        composition().zetav
                        (
                            speciei,
                            pCells[celli],
                            TvCells[celli]
                        );
                }
                else if (downgradeSingleTv)
                {
                    TvCells[celli] = TvCellsMix[celli];
                    hvelCells[celli] =
                        composition().HEvel
                        (
                            speciei,
                            pCells[celli],
                            TvCells[celli]
                        );
                    zetavCells[celli] =
                        composition().zetav
                        (
                            speciei,
                            pCells[celli],
                            TvCells[celli]
                        );
                }
                else if (composition().vibTempAssociativity(speciei) == -1)
                {
                    if
                    (
                        composition().noVibrationalTemp(speciei) == 1
                     or downgradeSingleVibMode
                    )
                    {
                        if (TtCells[celli] > vibrationalCutOffTemp)
                        {
                            if (YCells[celli] > miniYforSolvingEvEqn)
                            {
                                // it means that evEqn IS solved because Y_m
                                // large enough
                                
                                // TODO analyse and revise if necessary
                                if (sign(hvelCells[celli]) == -1) 
                                {
                                    //Info << "hvelCells[" << celli << "]" << tab << hvelCells[celli] << endl;
                                    //Info << "B1: TvCells[" << celli << "]" << tab << TvCells[celli] << endl;

                                    hvelCells[celli] =
                                        composition().HEvel
                                        (
                                            speciei,
                                            pCells[celli],
                                            vibrationalCutOffTemp
                                        );
                                    TvCells[celli] =
                                        TvelEvels
                                        (
                                            speciei,
                                            hvelCells[celli],
                                            pCells[celli],
                                            TvCells[celli]
                                        );
                                    //Info << "A1: TvCells[" << celli << "]" << tab << TvCells[celli] << endl;
                                }
                                else
                                {
                                    //Info << "B2: TvCells[" << celli << "]" << tab << TvCells[celli] << endl;
                                    TvCells[celli] =
                                        TvelEvels
                                        (
                                            speciei,
                                            hvelCells[celli],
                                            pCells[celli],
                                            TvCells[celli]
                                        );
                                    //Info << "A2: TvCells[" << celli << "]" << tab << TvCells[celli] << endl;
                                }
                            }
                            else
                            {
                                // TODO improve the treatment
                                TvCells[celli] = Tlow_; //TtCells[celli]; 
                                //Info << "A3: TvCells[" << celli << "]" << tab << TvCells[celli] << endl;
                            }
                        }
                        else
                        {
                            TvCells[celli] = TtCells[celli];
                            //Info << "A4: TvCells[" << celli << "]" << tab << TvCells[celli] << endl;
                        }

                        hvelCellsMix[celli] += YCells[celli]*hvelCells[celli];
                        zetavCells[celli] =
                            composition().zetav
                            (
                                speciei,
                                pCells[celli],
                                TvCells[celli]
                            );
                    }
                    else
                    {
                        // ABORTIVE WORK
                        /*zetavCells[celli] = 0.0;
                        TvCells[celli] = 0.0;
                        hvelCells[celli] = 0.0;

                        forAll(Tv_mode[speciei], mode)
                        {
                            scalarField& hvelmodeCells =
                                hevel_mode[speciei][mode].internalField();
                            scalarField& TvmodeCells =
                                Tv_mode[speciei][mode].internalField();
                            scalarField& zetavmodeCells =
                                zetav_mode[speciei][mode].internalField();
                            TvmodeCells[celli] =
                                composition().TvelHEvel_mode
                                (
                                    speciei,
                                    mode,
                                    hvelmodeCells[celli],
                                    pCells[celli],
                                    TvmodeCells[celli]
                                );
                            zetavmodeCells[celli] =
                                composition().zetav_mode
                                (
                                    speciei,
                                    mode,
                                    pCells[celli],
                                    TvmodeCells[celli]
                                );
                            zetavCells[celli] += zetavmodeCells[celli];
                            hvelCells[celli] += hvelmodeCells[celli];
                            TvCells[celli] +=
                                TvmodeCells[celli]*zetavmodeCells[celli];
                        }

                        TvCells[celli] /= zetavCells[celli];*/
                    }
                }
                else
                {
                    // If no molecule in the flow-field, then necessarily the
                    // 1-T solver is run using downgradeSingleTemperature
                    // because the electronic energy equation is not implemented
                    // and it cannot be substituted by a vibrational energy
                    // equation. Therefore, the following always works.
                    TvCells[celli] =
                        Tv[composition().vibTempAssociativity(speciei)]
                            .internalField()[celli];
                    hvelCells[celli] =
                        composition().HEvel
                        (
                            speciei,
                            pCells[celli],
                            TvCells[celli]
                        );
                    zetavCells[celli] =
                        composition().zetav
                        (
                            speciei,
                            pCells[celli],
                            TvCells[celli]
                        );
                    hvelCellsMix[celli] += YCells[celli]*hvelCells[celli];
                }

                hvCells[celli] =
                    composition().HEv(speciei, pCells[celli], TvCells[celli]);
                helCells[celli] =
                    composition().HEel(speciei, pCells[celli], TvCells[celli]);

                hvCellsMix[celli] += YCells[celli]*hvCells[celli];
                helCellsMix[celli] += YCells[celli]*helCells[celli];
                
                if (not downgradeSingleTv)
                {
                    TvCellsMix[celli] +=
                        XCells[celli]*zetavCells[celli]*TvCells[celli];
                }
                
                zetavCellsMix[celli] += XCells[celli]*zetavCells[celli];
                totZetavCellsMix[celli] += XCells[celli]*zetavCells[celli];
            }

            CvtrCells[celli] =
                composition().Cv_t(speciei, pCells[celli], TtCells[celli]);
            CvvelCells[celli] =
                composition().Cv_vel(speciei, pCells[celli], TvCells[celli]);
            CvtrCellsMix[celli] += YCells[celli]*CvtrCells[celli];
            CvvelCellsMix[celli] += YCells[celli]*CvvelCells[celli];
            CptrCellsMix[celli] += YCells[celli]*CvtrCells[celli];
            CpvelCellsMix[celli] += YCells[celli]*CvvelCells[celli];
            
            if (composition().isHeavySpecies(speciei))
            {
                CptrCellsMix[celli] += YCells[celli]*R;
                
                pPCells[celli] = YCells[celli]
                    /composition().psi(speciei, pCells[celli], TtCells[celli]);
            }
            else
            {
                CpvelCellsMix[celli] += YCells[celli]*R;
                
                pPCells[celli] = YCells[celli]
                    /composition().psi(speciei, pCells[celli], TvCells[celli]);
            }
            
            psiCellsMix[celli] += pPCells[celli];
            pPCells[celli] *= rhoCells[celli];
            
        }//end cells loop
    }//end species loop
    
    forAll(pCells, celli)
    {
        psiCellsMix[celli] = 1.0/psiCellsMix[celli];
        
        CvCellsMix[celli] = CvtrCellsMix[celli] + CvvelCellsMix[celli];
        CpCellsMix[celli] = CvCellsMix[celli] + RmixCells[celli];
        
        if (not downgradeSingleTv and totZetavCellsMix[celli] > SMALL)
        {
            TvCellsMix[celli] /= totZetavCellsMix[celli];
        }
    }

    //- Patch values calculation
    forAll(this->T().boundaryField(), patchi)
    {
        const fvPatchScalarField& pp = this->p().boundaryField()[patchi];
        const fvPatchScalarField& prho = this->rho_.boundaryField()[patchi];
        const fvPatchScalarField& pTt = this->T().boundaryField()[patchi];
        const fvPatchScalarField& pRmix = this->RMix_.boundaryField()[patchi];

//        fvPatchScalarField& phMix =
//            composition().e().boundaryFieldRef()[patchi];  // NEW VINCENT 2021/01/05
//        const fvPatchScalarField& phtMix = this->het().boundaryField()[patchi];  // NEW VINCENT 2021/01/05
        fvPatchScalarField& phvMix = this->hevMix_.boundaryFieldRef()[patchi];
        fvPatchScalarField& phelMix = this->heelMix_.boundaryFieldRef()[patchi];
        fvPatchScalarField& phvelMix = this->hevel().boundaryFieldRef()[patchi];
        fvPatchScalarField& pTvMix = this->Tv().boundaryFieldRef()[patchi];
        fvPatchScalarField& pzetavMix = this->zetav_.boundaryFieldRef()[patchi];
        fvPatchScalarField ptotZetavMix = pzetavMix;

        fvPatchScalarField& ppsiMix = this->psi_.boundaryFieldRef()[patchi];
        fvPatchScalarField& pCvtrMix =
            this->CvtrMix_.boundaryFieldRef()[patchi];
        fvPatchScalarField& pCvvelMix =
            this->CvvelMix_.boundaryFieldRef()[patchi];
        fvPatchScalarField& pCptrMix =
            this->CptrMix_.boundaryFieldRef()[patchi];
        fvPatchScalarField& pCpvelMix =
            this->CpvelMix_.boundaryFieldRef()[patchi];
        fvPatchScalarField& pCvMix = this->CvMix_.boundaryFieldRef()[patchi];
        fvPatchScalarField& pCpMix = this->CpMix_.boundaryFieldRef()[patchi];

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
            pCvtrMix[facei] = 0.0;
            pCvvelMix[facei] = 0.0;
            pCptrMix[facei] = 0.0;
            pCpvelMix[facei] = 0.0;
        }

        if (pTt.fixesValue())
        {
            // the temperature is fixed ... the energy is calculated at patches
            forAll(Y, speciei)
            {
                const fvPatchScalarField& pY =
                    Y[speciei].boundaryField()[patchi];
                const fvPatchScalarField& pX =
                    X[speciei].boundaryField()[patchi];
                fvPatchScalarField& pTv =
                    Tv[speciei].boundaryFieldRef()[patchi];
                fvPatchScalarField& phv =
                    hv[speciei].boundaryFieldRef()[patchi];
                fvPatchScalarField& phel =
                    hel[speciei].boundaryFieldRef()[patchi];
                fvPatchScalarField& phvel =
                    hvel[speciei].boundaryFieldRef()[patchi];
                fvPatchScalarField& pzetav =
                    zetav[speciei].boundaryFieldRef()[patchi];

                forAll(pTt, facei)
                {
                    if (pTv[facei] > SMALL)
                    {
                        if (downgradeSingleTemperature)
                        {
                            pTv[facei] = pTt[facei]; // NEW VINCENT 2021/01/05
                            phv[facei] =
                                composition().HEv
                                (
                                    speciei,
                                    pp[facei],
                                    pTv[facei]
                                );
                            phel[facei] =
                                composition().HEel
                                (
                                    speciei,
                                    pp[facei],
                                    pTv[facei]
                                );
                            phvel[facei] =
                                composition().HEvel
                                (
                                    speciei,
                                    pp[facei],
                                    pTv[facei]
                                );
                            phvelMix[facei] += pY[facei]*phvel[facei];
                            pzetav[facei] =
                                composition().zetav
                                (
                                    speciei,
                                    pp[facei],
                                    pTv[facei]
                                );
//                            pTvMix[facei] += pX[facei]*pzetav[facei]*pTv[facei]; // NEW VINCENT 2021/01/05
                        }
                        else if (downgradeSingleTv)
                        {
                            pTv[facei] = pTvMix[facei];
                            phv[facei] =
                                composition().HEv
                                (
                                    speciei,
                                    pp[facei],
                                    pTv[facei]
                                );
                            phel[facei] =
                                composition().HEel
                                (
                                    speciei,
                                    pp[facei],
                                    pTv[facei]
                                );
                            phvel[facei] =
                                composition().HEvel
                                (
                                    speciei,
                                    pp[facei],
                                    pTv[facei]
                                );
                            pzetav[facei] =
                                composition().zetav
                                (
                                    speciei,
                                    pp[facei],
                                    pTv[facei]
                                );
                        }
                        else
                        {
//                            phvel[facei] =
//                                composition().HEvel
//                                (
//                                    speciei,
//                                    pp[facei],
//                                    pTv[facei]
//                                );
                            phv[facei] =
                                composition().HEv
                                (
                                    speciei,
                                    pp[facei],
                                    pTv[facei]
                                );
                            phel[facei] =
                                composition().HEel
                                (
                                    speciei,
                                    pp[facei],
                                    pTv[facei]
                                );
                            phvelMix[facei] += pY[facei]*phvel[facei];
                            pzetav[facei] =
                                composition().zetav
                                (
                                    speciei,
                                    pp[facei],
                                    pTv[facei]
                                );
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
//                phMix[facei] = phtMix[facei] + phvelMix[facei]; // NEW VINCENT 2021/01/05

                if (downgradeSingleTemperature)
                {
                    pTvMix[facei] = pTt[facei];
                }
                else if (not downgradeSingleTv and ptotZetavMix[facei] > SMALL)
                {
                    pTvMix[facei] /= ptotZetavMix[facei];
                }
            }
        }
        else
        {
            // condition on the energy fields ... the temperatures are
            // calculated at patches
            forAll(pTt, facei)
            {
                forAll(Y, speciei)
                {
                    const fvPatchScalarField& pY =
                        Y[speciei].boundaryField()[patchi];
                    const fvPatchScalarField& pX =
                        X[speciei].boundaryField()[patchi];
                    const fvPatchScalarField& phv =
                        hv[speciei].boundaryField()[patchi];
                    const fvPatchScalarField& phel =
                        hel[speciei].boundaryField()[patchi];

                    fvPatchScalarField& phvel =
                        hvel[speciei].boundaryFieldRef()[patchi];
                    fvPatchScalarField& pTv =
                        Tv[speciei].boundaryFieldRef()[patchi];
                    fvPatchScalarField& pzetav =
                        zetav[speciei].boundaryFieldRef()[patchi];

                    if (pTv[facei] > SMALL)
                    {
                        if (downgradeSingleTemperature)
                        {
                            pTv[facei] = pTt[facei];
                            phvel[facei] =
                                composition().HEvel
                                (
                                    speciei,
                                    pp[facei],
                                    pTv[facei]
                                );
                            pzetav[facei] =
                                composition().zetav
                                (
                                    speciei,
                                    pp[facei],
                                    pTv[facei]
                                );
                            phvelMix[facei] += pY[facei]*phvel[facei];
//                            pTvMix[facei] += pX[facei]*pzetav[facei]*pTv[facei]; // NEW VINCENT 2021/01/05
                        }
                        else if (downgradeSingleTv)
                        {
                            pTv[facei] = pTvMix[facei];
                            phvel[facei] =
                                composition().HEvel
                                (
                                    speciei,
                                    pp[facei],
                                    pTv[facei]
                                );
                            pzetav[facei] =
                                composition().zetav
                                (
                                    speciei,
                                    pp[facei],
                                    pTv[facei]
                                );
                        }
                        else if (composition().vibTempAssociativity(speciei) == -1)
                        {
                            if
                            (
                                composition().noVibrationalTemp(speciei) == 1
                             or downgradeSingleVibMode
                            )
                            {
                                if (pTt[facei] > vibrationalCutOffTemp)
                                {
                                    if (pY[facei] > miniYforSolvingEvEqn)
                                    {
                                        // itmeans that evEqn IS solved because
                                        // Y_m large enough
                                        
                                        // NEW VINCENT 09/02/2017
                                        if (sign(phvel[facei]) == -1) 
                                        {
                                            phvel[facei] =
                                                composition().HEvel
                                                (
                                                    speciei,
                                                    pp[facei],
                                                    vibrationalCutOffTemp
                                                );
                                        }

                                        pTv[facei] =
                                            TvelEvels
                                            (
                                                speciei,
                                                phvel[facei],
                                                pp[facei],
                                                pTv[facei]
                                            );
                                    }
                                    else
                                    {
                                        // TODO improve the treatment
                                        pTv[facei] = TlowPatches_; //pTt[facei];
                                    }
                                }

                                pzetav[facei] =
                                    composition().zetav
                                    (
                                        speciei,
                                        pp[facei],
                                        pTv[facei]
                                    );
                                phvelMix[facei] += pY[facei]*phvel[facei];
                                pTvMix[facei] += pX[facei]
                                    *pzetav[facei]*pTv[facei];
                            }
                            else
                            {
                                // ABORTIVE WORK
                                /*pzetav[facei] = 0.0;
                                pTv[facei] = 0.0;

                                forAll(Tv_mode[speciei], mode)
                                {
                                    Tv_mode[speciei][mode].boundaryField()[patchi][facei] =
                                        composition().TvelHEvel_mode
                                        (
                                            speciei,
                                            mode,
                                            phvel[facei],
                                            pp[facei],
                                            Tv_mode[speciei][mode].boundaryField()[patchi][facei]
                                        );
                                    zetav_mode[speciei][mode].boundaryField()[patchi][facei] =
                                        composition().zetav_mode
                                        (
                                            speciei,
                                            mode,
                                            pp[facei],
                                            Tv_mode[speciei][mode]
                                                .boundaryField()[patchi][facei]
                                        );
                                    pzetav[facei] +=
                                        zetav_mode[speciei][mode]
                                            .boundaryField()[patchi][facei];
                                    pTv[facei] +=
                                        Tv_mode[speciei][mode]
                                            .internalField()[facei]
                                      * zetav_mode[speciei][mode]
                                            .boundaryField()[patchi][facei];
                                }

                                pTv[facei] /= pzetav[facei];*/
                            }
                        }
                        else
                        {
                            const fvPatchScalarField& pTvMol =
                                Tv[composition().vibTempAssociativity(speciei)]
                                    .boundaryField()[patchi];
                            pTv[facei] = pTvMol[facei];
                            pzetav[facei] =
                                composition().zetav
                                (
                                    speciei,
                                    pp[facei],
                                    pTv[facei]
                                );
                        }

                        phvMix[facei] += pY[facei]*phv[facei];
                        phelMix[facei] += pY[facei]*phel[facei];
                        pzetavMix[facei] += pX[facei]*pzetav[facei];
                        ptotZetavMix[facei] += pX[facei]*pzetav[facei];
                    }
                }//end species loop

                if (downgradeSingleTemperature)
                {
                    pTvMix[facei] = pTt[facei];
                }
                else if (not downgradeSingleTv and ptotZetavMix[facei] > SMALL)
                {
                    pTvMix[facei] /= ptotZetavMix[facei];
                }
            }//end faces loop
        }

        forAll(Y, speciei)
        {
            const scalar R = composition().R(speciei);
            const fvPatchScalarField& pY = Y[speciei].boundaryField()[patchi];
            const fvPatchScalarField& pTv = Tv[speciei].boundaryField()[patchi];
            
            fvPatchScalarField& ppP = pP[speciei].boundaryFieldRef()[patchi];
            fvPatchScalarField& pCvtr =
                Cvtr[speciei].boundaryFieldRef()[patchi];
            fvPatchScalarField& pCvvel =
                Cvvel[speciei].boundaryFieldRef()[patchi];

            forAll(pTt, facei)
            {
                pCvtr[facei] =
                    composition().Cv_t(speciei, pp[facei], pTt[facei]);
                pCvvel[facei] =
                    composition().Cv_vel(speciei, pp[facei], pTv[facei]);
                pCvtrMix[facei] += pY[facei]*pCvtr[facei];
                pCvvelMix[facei] += pY[facei]*pCvvel[facei];
                pCptrMix[facei] += pY[facei]*pCvtr[facei];
                pCpvelMix[facei] += pY[facei]*pCvvel[facei];

                if (composition().isHeavySpecies(speciei))
                {
                    pCptrMix[facei] += pY[facei]*R;
                    
                    ppP[facei] = pY[facei]
                        /composition().psi(speciei, pp[facei], pTt[facei]);
                }
                else
                {
                    pCpvelMix[facei] += pY[facei]*R;
                    
                    ppP[facei] = pY[facei]
                        /composition().psi(speciei, pp[facei], pTv[facei]);
                }
                
                ppsiMix[facei] += ppP[facei];
                ppP[facei] *= prho[facei];
                
            }//end faces loop
        }//end species loop
        
        ppsiMix = 1.0/ppsiMix;
        
        pCvMix = pCvtrMix + pCvvelMix;
        pCpMix = pCvMix + pRmix;
        
    }//end patches loop
    
    // Calculation of the electron pressure 
    if (composition().contains("e-"))
    {
        this->pe_ = composition().pP("e-");
    }
}


// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //


void Foam::rho2ReactionThermo::calculateLight()
{
    // TODO
}


// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //


void Foam::rho2ReactionThermo::calculateHEVK()
{
    // Tve is set according to the vibTempAssociativity table for particles
    // whose vibro-electronic energy equation is not solved

    //- Declarations
    const bool downgradeSingleTv = this->downgradeSingleTv();

    const scalarField& pCells = this->p().internalField();
    const scalarField& TvCellsMix = this->Tv().internalField();
    PtrList<Foam::volScalarField>& Tv = composition().Tv();
    PtrList<Foam::volScalarField>& hvel = composition().hevel();

    //- Cell values
    forAll(Tv, speciei)
    {
        scalarField& TvCells = Tv[speciei].primitiveFieldRef();
        scalarField& hvelCells = hvel[speciei].primitiveFieldRef();

        if (TvCells[0] > SMALL)
        {
            if (downgradeSingleTv)
            {
                forAll(pCells, celli)
                {
                    TvCells[celli] = TvCellsMix[celli];
                    hvelCells[celli] =
                        composition().HEvel
                        (
                            speciei,
                            pCells[celli],
                            TvCells[celli]
                        );
                }
            }
            else if (composition().vibTempAssociativity(speciei) != -1)
            {
                forAll(pCells, celli)
                {
                    TvCells[celli] =
                        Tv[composition().vibTempAssociativity(speciei)]
                            .internalField()[celli];
                    hvelCells[celli] =
                        composition().HEvel
                        (
                            speciei,
                            pCells[celli],
                            TvCells[celli]
                        );
                }
            }
        }

        //- Patch values
        forAll(this->p().boundaryField(), patchi)
        {
            fvPatchScalarField& pTv = Tv[speciei].boundaryFieldRef()[patchi];
            const fvPatchScalarField& pTvMix =
                this->Tv().boundaryField()[patchi];

            if (pTv.size() != 0)
            {
                if (pTv[0] > SMALL)
                {
                    const fvPatchScalarField& pp =
                        this->p().boundaryField()[patchi];
                    fvPatchScalarField& phvel =
                        hvel[speciei].boundaryFieldRef()[patchi];

                    if (downgradeSingleTv)
                    {
                        forAll(pTv, facei)
                        {
                            pTv[facei] = pTvMix[facei];
                            phvel[facei] =
                                composition().HEvel
                                (
                                    speciei,
                                    pp[facei],
                                    pTv[facei]
                                );
                        }
                    }
                    else if (composition().vibTempAssociativity(speciei) != -1)
                    {
                        forAll(pTv, facei)
                        {
                            pTv[facei] =
                                Tv[composition().vibTempAssociativity(speciei)]
                                    .boundaryField()[patchi][facei];
                            phvel[facei] =
                                composition().HEvel
                                (
                                    speciei,
                                    pp[facei],
                                    pTv[facei]
                                );
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
    
    const volScalarField::Boundary wallPatches = this->T().boundaryField();

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
    if (this->downgradeSingleTemperature())
    {
        het_ = new volScalarField
        (
            IOobject
            (
                this->phasePropertyName("het"),
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimEnergy/dimMass
        );
    }
    else
    {
        het_ = new volScalarField
        (
            IOobject
            (
                this->phasePropertyName("het"),
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimEnergy/dimMass,
            this->het2BoundaryTypes(),
            this->het2BoundaryBaseTypes()
        );
    }
        
    if (this->downgradeSingleTv())
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


void Foam::rho2ReactionThermo::correctVolFractions()
{
    if (debug)
    {
        Info<< "entering rho2ReactionThermo::correctVolFractions()" << endl;
    }

    correctVolChemFractions();

    if (debug)
    {
        Info<< "exiting rho2ReactionThermo::correctVolFractions()" << endl;
    }
}


void Foam::rho2ReactionThermo::correctBdrFractions()
{
    if (debug)
    {
        Info<< "entering rho2ReactionThermo::correctBdrFractions()" << endl;
    }

    correctBdrChemFractions();

    if (debug)
    {
        Info<< "exiting rho2ReactionThermo::correctBdrFractions()" << endl;
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


Foam::word Foam::rho2ReactionThermo::partialThermoName()
{
    return partialThermoName_;
}


// ************************************************************************* //
