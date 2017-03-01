/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
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

const Foam::scalar Foam::rho2ReactionThermo::vibrationalCutOffTemp = 200;

const Foam::scalar Foam::rho2ReactionThermo::miniYforSolvingEvEqn = 1.0e-4;

bool Foam::rho2ReactionThermo::temperatureFieldOutOfRange = false;

bool Foam::rho2ReactionThermo::hasCrashedButRecovered = false;

Foam::label Foam::rho2ReactionThermo::noCellsWithTemperatureFieldOutOfRange = 0;

Foam::FixedList<Foam::scalar, 2> Foam::rho2ReactionThermo::minMaxTemperatureFieldOutOfRange;

// * * * * * * * * * * * * * private Member Functions  * * * * * * * * * * * //

Foam::word Foam::rho2ReactionThermo::transportToTypedef(const word transportModel)
{
    word typeDefName = word::null;
    
    if(transportModel == "constant")
    {
        typeDefName = "demConstGasEThermoPhysicsH2TGD";
    }
    else if(transportModel == "BlottnerEucken")
    {
        typeDefName = "demBEGasEThermoPhysicsH2TGD";
    }
    else if(transportModel == "powerLawEucken")
    {
        typeDefName = "demPLEGasEThermoPhysicsH2TGD";
    }
    else if(transportModel == "SutherlandEucken")
    {
        typeDefName = "demGasEThermoPhysicsH2TGD";
    }
    else if(transportModel == "CEA")
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
    
    scalarField sumX(0*X[0].internalField()), sumnD(0*nD[0].internalField()), 
        sumpP(0*pP[0].internalField()), sumpD(0*pD[0].internalField());
    
    forAll(Y, speciei)
    {           
        const scalarField& YCells = Y[speciei].internalField();
        
        scalarField& XCells = X[speciei].internalField();
        scalarField& nDCells = nD[speciei].internalField();
        scalarField& pPCells = pP[speciei].internalField();
        scalarField& pDCells = pD[speciei].internalField();
        
        if(speciei < Y.size() - 1) // NEW VINCENT 25/04/2016
        // This loop ensures that the sum of the chemical quantities are bounded 
        // (needs to be < Y.size()-1 to be activated)
        { 
            forAll(pCells, celli)
            { 
                XCells[celli] = composition().molarFraction(speciei, YCells[celli], celli);
                pPCells[celli] = composition().partialPressure(XCells[celli], pCells[celli]);
                nDCells[celli] = composition().numberDensity(speciei, YCells[celli], rhoCells[celli]);
                pDCells[celli] = composition().partialDensity(YCells[celli], rhoCells[celli]);
                
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
                XCells[celli] = max(1 - sumX[celli], 0); 
                nDCells[celli] = composition().numberDensity(speciei, YCells[celli], rhoCells[celli]);
                pPCells[celli] = max(pCells[celli] - sumpP[celli], 0);
                pDCells[celli] = max(rhoCells[celli] - sumpD[celli], 0);
            }//end cells loop
        }  
    }//end species loop 
    
    
    forAll(this->p_.boundaryField(), patchi)
    {
        const fvPatchScalarField& pp  = this->p_.boundaryField()[patchi];
        const fvPatchScalarField& prho = this->rho_.boundaryField()[patchi];
        
        scalarField sumX(0*X[0].boundaryField()[patchi]), sumnD(0*nD[0].boundaryField()[patchi]), 
            sumpP(0*pP[0].boundaryField()[patchi]), sumpD(0*pD[0].boundaryField()[patchi]);
        
        forAll(Y, speciei)
        {
            const fvPatchScalarField& pY = Y[speciei].boundaryField()[patchi];
            
            fvPatchScalarField& pX = X[speciei].boundaryField()[patchi];
            fvPatchScalarField& pnD = nD[speciei].boundaryField()[patchi];
            fvPatchScalarField& ppP = pP[speciei].boundaryField()[patchi];
            fvPatchScalarField& ppD = pD[speciei].boundaryField()[patchi]; 
            
            if(speciei < Y.size() - 1) // NEW VINCENT 25/04/2016
            {
                forAll(pp, facei)
                {
                    pX[facei]  = composition().molarFraction(speciei, pY[facei], patchi, facei);
                    pnD[facei] = composition().numberDensity(speciei, pY[facei], prho[facei]);
                    ppP[facei] = composition().partialPressure(pX[facei], pp[facei]);
                    ppD[facei] = composition().partialDensity(pY[facei], prho[facei]);
                    
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
                    pX[facei] = max(1 - sumX[facei], 0);
                    pnD[facei] = composition().numberDensity(speciei, pY[facei], prho[facei]);
                    ppP[facei] = max(pp[facei] - sumpP[facei], 0);
                    ppD[facei] = max(prho[facei] - sumpD[facei], 0);
                }//end faces loop 
            }      
        }//end species loop 
    }//end patches loop         
}


void Foam::rho2ReactionThermo::correctOverallTemperature()
{
    //- Declarations
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
    
    scalarField& zetarCellsMix = this->zetar_.internalField();
    scalarField& zetaelCellsMix = this->zetael_.internalField();
    scalarField& numTCellsMix = numTMix.internalField();
    scalarField& denTCellsMix = denTMix.internalField();
    scalarField& TCells = this->T_.internalField();
    
    //- Initialisations
    zetarCellsMix = 0.0;
    zetaelCellsMix = 0.0;
    numTCellsMix = 0.0;
    denTCellsMix = 0.0;
    
    forAll(this->p_.boundaryField(), patchi)
    {
        fvPatchScalarField& pzetarMix = this->zetar_.boundaryField()[patchi];
        fvPatchScalarField& pzetaelMix = this->zetael_.boundaryField()[patchi];
        fvPatchScalarField& pnumTMix = numTMix.boundaryField()[patchi];
        fvPatchScalarField& pdenTMix = denTMix.boundaryField()[patchi];
            
        pzetarMix = 0.0;
        pzetaelMix = 0.0;
        pnumTMix = 0.0;
        pdenTMix = 0.0;
    }
        
    //- Field calculations
    forAll(Y, speciei)
    {
        //- Cells values
        const scalarField& YCells = Y[speciei].internalField();
        const scalarField& XCells = X[speciei].internalField();
        const scalarField& TvCells = Tv[speciei].internalField();
        const scalarField& zetavCells = zetav[speciei].internalField();
        scalarField& zetarCells = zetar[speciei].internalField();
        scalarField& zetaelCells = zetael[speciei].internalField();
        
        forAll(pCells, celli)
        {           
            zetarCells[celli] = composition().zetar(speciei, pCells[celli], TtCells[celli], TvCells[celli]);
            zetaelCells[celli] = composition().zetael(speciei, pCells[celli], TvCells[celli]);
            zetarCellsMix[celli] += XCells[celli]*zetarCells[celli];
            zetaelCellsMix[celli] += XCells[celli]*zetaelCells[celli];
            
            if(composition().species()[speciei] == "e-")
            {
                numTCellsMix[celli] += 3.0*TvCells[celli]*YCells[celli];
                denTCellsMix[celli] += 3.0*YCells[celli];
            }
            else
            {
                numTCellsMix[celli] += ((3.0 + zetarCells[celli])*TtCells[celli] + (zetavCells[celli] + zetaelCells[celli])*TvCells[celli])*YCells[celli];
                denTCellsMix[celli] += (3.0 + zetarCells[celli] + zetavCells[celli] + zetaelCells[celli])*YCells[celli]; 
            }
        }//end cells loop
    
        //- Patch values calculations
        forAll(this->p_.boundaryField(), patchi)
        {
            const fvPatchScalarField& pp = this->p_.boundaryField()[patchi];
            const fvPatchScalarField& pTt = this->Tt_.boundaryField()[patchi];
            const fvPatchScalarField& pY = Y[speciei].boundaryField()[patchi];
            const fvPatchScalarField& pX = X[speciei].boundaryField()[patchi];
            const fvPatchScalarField& pTv = Tv[speciei].boundaryField()[patchi];
            const fvPatchScalarField& pzetav = zetav[speciei].boundaryField()[patchi];
            
            fvPatchScalarField& pzetar = zetar[speciei].boundaryField()[patchi];
            fvPatchScalarField& pzetael = zetael[speciei].boundaryField()[patchi];
            
            fvPatchScalarField& pzetarMix = this->zetar_.boundaryField()[patchi];
            fvPatchScalarField& pzetaelMix = this->zetael_.boundaryField()[patchi];
            fvPatchScalarField& pnumTMix = numTMix.boundaryField()[patchi];
            fvPatchScalarField& pdenTMix = denTMix.boundaryField()[patchi];
            
            forAll(pTt, facei) 
            {
                pzetar[facei] = composition().zetar(speciei, pp[facei], pTt[facei], pTv[facei]);
                pzetael[facei] = composition().zetael(speciei, pp[facei], pTv[facei]);
                pzetarMix[facei] += pX[facei]*pzetar[facei];
                pzetaelMix[facei] += pX[facei]*pzetael[facei];
                
                if(composition().species()[speciei] == "e-")
                {
                    pnumTMix[facei] += 3.0*pTv[facei]*pY[facei];
                    pdenTMix[facei] += 3.0*pY[facei];
                }
                else
                {
                    pnumTMix[facei] += ((3.0 + pzetar[facei])*pTt[facei] + (pzetav[facei] + pzetael[facei])*pTv[facei])*pY[facei];
                    pdenTMix[facei] += (3.0 + pzetar[facei] + pzetav[facei] + pzetael[facei])*pY[facei];
                }
            }//end faces loop
        }//end patches loop
    }//end species loop
    
    forAll(pCells, celli)
    {
        TCells[celli] = numTCellsMix[celli]/denTCellsMix[celli];
    }
    
    forAll(this->p_.boundaryField(), patchi)
    {
        const fvPatchScalarField& pnumTMix = numTMix.boundaryField()[patchi];
        const fvPatchScalarField& pdenTMix = denTMix.boundaryField()[patchi];
        fvPatchScalarField& pT = this->T_.boundaryField()[patchi];
        
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
    
    //- Declarations
    const bool downgradeSingleTemperature = this->downgradeSingleTemperature();
    const bool downgradeSingleTv = this->downgradeSingleTv();
    const PtrList<Foam::volScalarField>& Y = composition().Y();
    const PtrList<Foam::volScalarField>& X = composition().X();
    const scalarField& pCells = this->p_.internalField();
    const scalarField& TtCells = this->Tt_.internalField();
    
    PtrList<Foam::volScalarField>& hvel = composition().hevel();
    PtrList<Foam::volScalarField>& Tv = composition().Tv();
    
    const scalarField& TvCellsMix = this->Tv_.internalField();
    scalarField& htCellsMix = this->het_.internalField();
    scalarField& hvelCellsMix = this->hevel().internalField();
    scalarField& hCellsMix = composition().e().internalField();
    scalarField& psiCellsMix = this->psi_.internalField();
    scalarField& rhoCellsMix = this->rho_.internalField();
    
    PtrList<scalar> YList(Y.size()); // NEW VINCENT 14/02/2017

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
            
            scalarField& hvelCells = hvel[speciei].internalField();
            scalarField& TvCells = Tv[speciei].internalField();
            
            if (TvCells[celli] != 0.0)
            {
                if(downgradeSingleTemperature)
                {
                    TvCells[celli] = TtCells[celli];
                }
                else if(downgradeSingleTv)
                {
                    TvCells[celli] = TvCellsMix[celli];
                }
                
                hvelCells[celli] = composition().HEvel(speciei, pCells[celli], TvCells[celli]);
                hvelCellsMix[celli] += YCells[celli]*hvelCells[celli];
                
                hCellsMix[celli] += YCells[celli]*hvelCells[celli];
            }
            
            if(composition().particleType(speciei) > 0) 
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
        fvPatchScalarField& pTvMix = this->Tv_.boundaryField()[patchi];
        
        fvPatchScalarField& phtMix = this->het_.boundaryField()[patchi];
        fvPatchScalarField& phvelMix = this->hevel().boundaryField()[patchi];
        fvPatchScalarField& phMix = composition().e().boundaryField()[patchi];
        
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
                    fvPatchScalarField& pTv = Tv[speciei].boundaryField()[patchi];
                    
                    fvPatchScalarField& phvel = hvel[speciei].boundaryField()[patchi];
                    
                    phtMix[facei] += pY[facei]*composition().HEt(speciei, pp[facei], pTt[facei]);
                    
                    if (pTv[facei] != 0.0)
                    {
                        if(downgradeSingleTv)
                        {
                            pTv[facei] = pTvMix[facei];
                        }
                        
                        phvel[facei] = composition().HEvel(speciei, pp[facei], pTv[facei]);    
                        phvelMix[facei] += pY[facei]*phvel[facei]; 
                          
                        phMix[facei] += phvelMix[facei];
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
                    fvPatchScalarField& phvel = hvel[speciei].boundaryField()[patchi];
                    fvPatchScalarField& pTv = Tv[speciei].boundaryField()[patchi];
                    
                    phtMix[facei] += pY[facei]*composition().HEt(speciei, pp[facei], pTt[facei]);
                    
                    if(pTv[facei] != 0.0)
                    {
                        if(downgradeSingleTemperature)
                        {
                            phvel[facei] = composition().HEvel(speciei, pp[facei], pTt[facei]); // NEW VINCENT 14/02/2017
                            pTv[facei] = pTt[facei];
                        }
                        else if(downgradeSingleTv)
                        {
                            phvel[facei] = composition().HEvel(speciei, pp[facei], pTvMix[facei]); // NEW VINCENT 14/02/2017
                            pTv[facei] = pTvMix[facei];
                        }
                        else if(composition().vibTempAssociativity(speciei) == -1)
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
                }
                
                phMix[facei] += phtMix[facei];
            }
        }
        
        fvPatchScalarField& ppsiMix = this->psi_.boundaryField()[patchi];
        fvPatchScalarField& prhoMix = this->rho_.boundaryField()[patchi];
        
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
                if(composition().particleType(speciei) > 0) 
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
}


// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -//


void Foam::rho2ReactionThermo::initialise()
{
    correctChemFractions();
    
    //- Declarations
    const bool downgradeSingleTemperature = this->downgradeSingleTemperature();
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
    
    scalarField& htCellsMix = this->het_.internalField(); // NEW VINCENT
    scalarField& hvCellsMix = this->hevMix_.internalField();
    scalarField& helCellsMix = this->heelMix_.internalField();
    scalarField& hvelCellsMix = this->hevel().internalField();
    scalarField& hCellsMix = composition().e().internalField(); // NEW VINCENT 23/02/2016
    scalarField& TvCellsMix = this->Tv_.internalField();
    scalarField& zetavCellsMix = this->zetav_.internalField();
    scalarField totZetavCellsMix = zetavCellsMix; // NEW VINCENT 25/04/2016
    scalarField& psiCellsMix = this->psi_.internalField(); // NEW VINCENT 11/04/2016
    scalarField& rhoCellsMix = this->rho_.internalField(); // NEW VINCENT 11/04/2016

    //- Cells values
    forAll(pCells, celli)
    { 
        htCellsMix[celli] = 0.0; // NEW VINCENT
        hvCellsMix[celli] = 0.0;
        helCellsMix[celli] = 0.0;
        hvelCellsMix[celli] = 0.0;
        hCellsMix[celli] = 0.0; // NEW VINCENT
        TvCellsMix[celli] = 0.0;
        zetavCellsMix[celli] = 0.0;
        totZetavCellsMix[celli] = 0.0; // NEW VINCENT 25/04/2016
        
        psiCellsMix[celli] = 0.0; // NEW VINCENT 11/04/2016
        rhoCellsMix[celli] = 0.0; // NEW VINCENT 11/04/2016
        
        forAll(Y, speciei) // NEW VINCENT 05/03/2016
        {           
            const scalarField& YCells = Y[speciei].internalField();
            htCellsMix[celli] += YCells[celli]*composition().HEt(speciei, pCells[celli], TtCells[celli]);
        }
        
        forAll(Y, speciei)
        {           
            const scalarField& YCells = Y[speciei].internalField();
            const scalarField& XCells = X[speciei].internalField();
            
            scalarField& hvCells = hv[speciei].internalField();
            scalarField& helCells = hel[speciei].internalField();
            scalarField& hvelCells = hvel[speciei].internalField();
            scalarField& TvCells = Tv[speciei].internalField();
            scalarField& zetavCells = zetav[speciei].internalField();
            
            if (TvCells[celli] != 0.0)
            {
                hvCells[celli] = composition().HEv(speciei, pCells[celli], TvCells[celli]);
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
                
                helCells[celli] = composition().HEel(speciei, pCells[celli], TvCells[celli]);
                zetavCells[celli] = composition().zetav(speciei, pCells[celli], TvCells[celli]);     
                
                hvCellsMix[celli] += YCells[celli]*hvCells[celli];
                helCellsMix[celli] += YCells[celli]*helCells[celli];
                hvelCellsMix[celli] += YCells[celli]*hvelCells[celli];
                hCellsMix[celli] += YCells[celli]*hvelCells[celli];
                TvCellsMix[celli] += YCells[celli]*zetavCells[celli]*TvCells[celli]; // NEW VINCENT 25/04/2016
                zetavCellsMix[celli] += XCells[celli]*zetavCells[celli];
                totZetavCellsMix[celli] += YCells[celli]*zetavCells[celli]; // NEW VINCENT 25/04/2016
            }
            
            if(composition().particleType(speciei) > 0) 
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
        
        if(totZetavCellsMix[celli] != 0.0)
        {
            TvCellsMix[celli] /= totZetavCellsMix[celli]; // NEW VINCENT 25/04/2016
        }
        
        hCellsMix[celli] += htCellsMix[celli]; // NEW VINCENT
    }
    
    //- patch values calculations
    forAll(this->Tt_.boundaryField(), patchi)
    {
        const fvPatchScalarField& pp = this->p_.boundaryField()[patchi];
        const fvPatchScalarField& pTt = this->Tt_.boundaryField()[patchi];

        fvPatchScalarField& phtMix = this->het_.boundaryField()[patchi]; // NEW VINCENT
        fvPatchScalarField& phvMix = this->hevMix_.boundaryField()[patchi];
        fvPatchScalarField& phelMix = this->heelMix_.boundaryField()[patchi];
        fvPatchScalarField& phvelMix = this->hevel().boundaryField()[patchi];
        fvPatchScalarField& phMix = composition().e().boundaryField()[patchi]; // NEW VINCENT 23/02/2016
        fvPatchScalarField& pTvMix = this->Tv_.boundaryField()[patchi];
        
        if(pTt.fixesValue())
        // the temperature is fixed ... the energy is calculated at patches
        {
            forAll(pTt, facei)
            {  
                phtMix[facei] = 0.0; // NEW VINCENT
                phvMix[facei] = 0.0;
                phelMix[facei] = 0.0;
                phvelMix[facei] = 0.0;
                phMix[facei] = 0.0; // NEW VINCENT
                pTvMix[facei] = 0.0;
                
                forAll(Y, speciei)
                {
                    const fvPatchScalarField& pY = Y[speciei].boundaryField()[patchi];
                    const fvPatchScalarField& pTv = Tv[speciei].boundaryField()[patchi];
                    
                    fvPatchScalarField& phv = hv[speciei].boundaryField()[patchi];
                    fvPatchScalarField& phel = hel[speciei].boundaryField()[patchi];
                    fvPatchScalarField& phvel = hvel[speciei].boundaryField()[patchi];
                    
                    phtMix[facei] += pY[facei]*composition().HEt(speciei, pp[facei], pTt[facei]); // NEW VINCENT 25/01/16
                    
                    if (pTv[facei] != 0.0)
                    {
                        phv[facei] = composition().HEv(speciei, pp[facei], pTv[facei]);               
                        phvel[facei] = composition().HEvel(speciei, pp[facei], pTv[facei]);    
                        phel[facei] = composition().HEel(speciei, pp[facei], pTv[facei]);
                        
                        phvMix[facei] += pY[facei]*phv[facei];
                        phelMix[facei] += pY[facei]*phel[facei];  
                        phvelMix[facei] += pY[facei]*phvel[facei];   
                        phMix[facei] += pY[facei]*phvel[facei]; // NEW VINCENT 25/01/16
                        pTvMix[facei] += pY[facei]*pTv[facei];
                    }
                }
                
                phMix[facei] += phtMix[facei]; // NEW VINCENT
            }
        }
        else 
        // condition on the energy fields ... the temperatures are calculated at patches
        {
            forAll(pTt, facei)
            {  
                phtMix[facei] = 0.0; // NEW VINCENT
                phvMix[facei] = 0.0;
                phelMix[facei] = 0.0;
                phvelMix[facei] = 0.0;
                phMix[facei] = 0.0; // NEW VINCENT
                pTvMix[facei] = 0.0;
                
                forAll(Y, speciei)
                {
                    const fvPatchScalarField& pY = Y[speciei].boundaryField()[patchi];
                    
                    fvPatchScalarField& phv = hv[speciei].boundaryField()[patchi];
                    fvPatchScalarField& phel = hel[speciei].boundaryField()[patchi];
                    fvPatchScalarField& phvel = hvel[speciei].boundaryField()[patchi];
                    fvPatchScalarField& pTv = Tv[speciei].boundaryField()[patchi];
                    
                    phtMix[facei] += pY[facei]*composition().HEt(speciei, pp[facei], pTt[facei]); // NEW VINCENT 25/01/16
                    
                    if (pTv[facei] != 0.0)
                    {
                        phv[facei] = composition().HEv(speciei, pp[facei], pTv[facei]);               
                        phvel[facei] = composition().HEvel(speciei, pp[facei], pTv[facei]);    
                        phel[facei] = composition().HEel(speciei, pp[facei], pTv[facei]);
                        
                        if (downgradeSingleTemperature) // NEW VINCENT 25/01/16
                        {
                            pTv[facei] = pTt[facei];
                        }
                        else if (composition().particleType(speciei) == 2) // molecules
                        {
                            pTv[facei] = TvelEvels(speciei, phvel[facei], pp[facei], pTv[facei]);
                            //composition().TvelHEvel(speciei, phvel[facei], pp[facei], pTv[facei]);
                        }
                        else // ions, electrons (or atoms if Eel on)
                        {
                            const fvPatchScalarField& pTvMol = Tv[composition().vibTempAssociativity(speciei)].boundaryField()[patchi];
                            pTv[facei] = pTvMol[facei];
                        }
                        
                        phvMix[facei] += pY[facei]*phv[facei];
                        phelMix[facei] += pY[facei]*phel[facei];
                        phvelMix[facei] += pY[facei]*phvel[facei];
                        phMix[facei] += pY[facei]*phvel[facei]; // NEW VINCENT 25/01/16
                        pTvMix[facei] += pY[facei]*pTv[facei];
                    }
                    else
                    {
                        phv[facei] = 0.0; 
                        phel[facei] = 0.0;
                        phvel[facei] = 0.0;           
                    }
                }
                
                phMix[facei] += phtMix[facei]; // NEW VINCENT
            }
        }
        
        fvPatchScalarField& ppsiMix = this->psi_.boundaryField()[patchi];
        fvPatchScalarField& prhoMix = this->rho_.boundaryField()[patchi];
        
        forAll(pTt, facei) 
        {
            ppsiMix[facei] = 0.0;
            prhoMix[facei] = 0.0;
            
            forAll(Y, speciei)
            {
                const fvPatchScalarField& pX = X[speciei].boundaryField()[patchi];
                const fvPatchScalarField& pTv = Tv[speciei].boundaryField()[patchi];
                
                if(composition().particleType(speciei) > 0) 
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
    PtrList<volScalarField>& fixedY,
    volScalarField& fixedTt,
    PtrList<volScalarField>& fixedTv,
    volScalarField& fixedRho,
    word& regionName
)
{
    label zoneID = this->Tt_.mesh().cellZones().findZoneID(regionName);

    if(zoneID == -1)
    {
        FatalErrorIn("fixedFields")
            << "Cannot find region: " << regionName << nl << "in: "
            << this->Tt_.mesh().time().constant()/"polyMesh/cellZones"
            << exit(FatalError);
    }
    
    const cellZone& zone = this->Tt_.mesh().cellZones()[zoneID];

    if (zone.size())
    {
        PtrList<Foam::volScalarField>& Y = composition().Y();
        
        scalarField& TtCells = this->Tt_.internalField();
        scalarField& fixedTtCells = fixedTt.internalField();
        PtrList<Foam::volScalarField>& Tv = composition().Tv();
        scalarField& rhoCells = this->rho_.internalField();
        scalarField& fixedRhoCells = fixedRho.internalField();
        
        scalarField& pCells = this->p_.internalField();
        scalarField& psiCells = this->psi_.internalField();
        scalarField& htCells = this->het_.internalField();
        PtrList<Foam::volScalarField>& hvel = composition().hevel();
        scalarField& hCells = composition().e().internalField();
        
        forAll(zone, idx)
        {
            const label celli = zone[idx];
            htCells[celli] = 0.0;
            hCells[celli] = 0.0;
            psiCells[celli] = 0.0;
            
            TtCells[celli] = fixedTtCells[celli];
            rhoCells[celli] = fixedRhoCells[celli];
        }
            
        forAll(Y, speciei)
        {
            scalarField& YCells = Y[speciei].internalField();
            scalarField& fixedYCells = fixedY[speciei].internalField();
            scalarField& TvCells = Tv[speciei].internalField();
            scalarField& fixedTvCells = fixedTv[speciei].internalField();
            scalarField& hvelCells = hvel[speciei].internalField();
                
            forAll(zone, idx)
            {
                const label celli = zone[idx];
                YCells[celli] = fixedYCells[celli];
                TvCells[celli] = fixedTvCells[celli];
                
                htCells[celli] += YCells[celli]*composition().HEt(speciei, pCells[celli], TtCells[celli]);
                if (TvCells[celli] != 0.0)
                {
                    hvelCells[celli] = composition().HEvel(speciei, pCells[celli], TvCells[celli]);
                }
                
                if(composition().particleType(speciei) > 0) 
                {
                    psiCells[celli] += composition().molarFraction(speciei, YCells[celli], celli)*composition().psi(speciei, pCells[celli], TtCells[celli]);
                }
                else
                {
                    psiCells[celli] += composition().molarFraction(speciei, YCells[celli], celli)*composition().psi(speciei, pCells[celli], TvCells[celli]);
                }
                
                hCells[celli] += hvelCells[celli];
            }
        }
        
        forAll(zone, idx)
        {
            const label celli = zone[idx];
            hCells[celli] += htCells[celli];
            pCells[celli] = rhoCells[celli]/psiCells[celli];
        }
        

        
       //TODO correct boundaries
    }
    
    correctChemFractions();
}


// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -//


void Foam::rho2ReactionThermo::calculate()
{
    //- Declarations
    const bool downgradeSingleTemperature = this->downgradeSingleTemperature();
    const bool downgradeSingleVibMode = this->downgradeSingleVibMode();
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
    
    scalarField& TtCells = this->Tt_.internalField();
    scalarField& TvCellsMix = this->Tv_.internalField();
    scalarField& hvCellsMix = this->hevMix_.internalField();
    scalarField& helCellsMix = this->heelMix_.internalField();
    scalarField& hvelCellsMix = this->hevel().internalField();
    scalarField& psiCellsMix = this->psi_.internalField();
    scalarField& zetavCellsMix = this->zetav_.internalField();
    scalarField totZetavCellsMix = zetavCellsMix;
    
    
    //- Initialisations
    TvCellsMix = 0.0;
    hvCellsMix = 0.0;
    helCellsMix = 0.0;
    if(not downgradeSingleTv) hvelCellsMix = 0.0;
    zetavCellsMix = 0.0;
    totZetavCellsMix = 0.0;
    psiCellsMix = 0.0;
    
       
    //- Overall temperature calculation (single-temperature model)
    //  or trans-rotational temperature calculation (two-temperature model)
    // NEW VINCENT 17/06/2016 *************************************************
    PtrList<scalar> YList (Y.size());

    forAll(pCells, celli)
    {
        forAll(YList, speciei)
        {
            YList.set(speciei, new scalar(Y[speciei].internalField()[celli]));
        } 

        if(downgradeSingleTemperature)
        {
            TtCells[celli] = TEs
            (
                hCellsMix[celli],
                pCells[celli],
                TtCells[celli],
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
    } 

    forAll(this->Tt_.boundaryField(), patchi)
    {
        const volScalarField::GeometricBoundaryField wallPatches = this->Tt_.boundaryField(); // NEW VINCENT 23/08/2016
        const fvPatchScalarField& pp = this->p_.boundaryField()[patchi];
        const fvPatchScalarField& phMix = composition().e().boundaryField()[patchi];
        const scalarField& phtMix = this->het_.boundaryField()[patchi];
        const scalarField& phvMix = this->hevel().boundaryField()[patchi];
        fvPatchScalarField& pTt = this->Tt_.boundaryField()[patchi];
        fvPatchScalarField& pTvMix = this->Tv_.boundaryField()[patchi];

        if (not pTt.fixesValue())
        // condition on the energy fields ... the temperature is calculated at patches
        {
            forAll(pTt, facei)
            {
                forAll(YList, speciei)
                {
                    YList[speciei] = Y[speciei].boundaryField()[patchi][facei];
                } 
            
                if(downgradeSingleTemperature)
                {
                    pTt[facei] = TEs
                    (
                        phMix[facei],
                        pp[facei],
                        pTt[facei],
                        YList
                    );
                }
                else if(downgradeSingleTv) // NEW VINCENT 16/08/2016
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
                        phvMix[facei],
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
            }
        }
        
        if(isA<wallFvPatch>(wallPatches[patchi].patch())) // NEW VINCENT 22/08/2016
        // to increase the stability of the 2T model in near-wall high Kn-number regions
        {
            forAll(pTt, facei)
            {
                if(pTt[facei] > ThighPatches_) pTt[facei] = ThighPatches_;
                else if(pTt[facei] < TlowPatches_) pTt[facei] = TlowPatches_;
            }
        }
    }
    // END NEW VINCENT 17/06/2016 *********************************************
    

    //- Cells values
    forAll(Y, speciei)
    {  
        const scalarField& YCells = Y[speciei].internalField();
        const scalarField& XCells = X[speciei].internalField();
        
        scalarField& TvCells = Tv[speciei].internalField();
        scalarField& hvCells = hv[speciei].internalField();
        scalarField& helCells = hel[speciei].internalField();
        scalarField& hvelCells = hvel[speciei].internalField();
        scalarField& zetavCells = zetav[speciei].internalField();
        
        forAll(pCells, celli)
        {
            if (TvCells[celli] != 0.0)
            {
                if(downgradeSingleTemperature)
                {
                    TvCells[celli] = TtCells[celli];
                    hvelCells[celli] = composition().HEvel(speciei, pCells[celli], TvCells[celli]);
                    zetavCells[celli] = composition().zetav(speciei, pCells[celli], TvCells[celli]);
                    hvelCellsMix[celli] += YCells[celli]*hvelCells[celli];
                }
                else if(downgradeSingleTv)
                {
                    // yet TODO
                    zetavCells[celli] = composition().zetav(speciei, pCells[celli], TvCells[celli]);
                }
                else if(composition().vibTempAssociativity(speciei) == -1)
                {
                    if(composition().noVibrationalTemp(speciei) == 1 or downgradeSingleVibMode)
                    {
                        if(TtCells[celli] > vibrationalCutOffTemp)
                        {
                            TvCells[celli] = TvelEvels(speciei, hvelCells[celli], pCells[celli], TvCells[celli]);
                            //composition().TvelHEvel(speciei, hvelCells[celli], pCells[celli], TvCells[celli]);
                        }
                        else
                        {
                            TvCells[celli] = TtCells[celli];
                        }
                        zetavCells[celli] = composition().zetav(speciei, pCells[celli], TvCells[celli]);
                        hvelCellsMix[celli] += YCells[celli]*hvelCells[celli];
                    }
                    else
                    {
                        /*zetavCells[celli] = 0.0; TODO ONGOING WORK 
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
                TvCellsMix[celli] += zetavCells[celli]*TvCells[celli];
                zetavCellsMix[celli] += XCells[celli]*zetavCells[celli];
                totZetavCellsMix[celli] += zetavCells[celli];
            }
            
            if(composition().particleType(speciei) > 0) 
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
        if(totZetavCellsMix[celli] != 0.0)
        {
            TvCellsMix[celli] /= totZetavCellsMix[celli]; // NEW VINCENT 25/04/2016
        }
    }

    //- Patch values
    forAll(this->Tt_.boundaryField(), patchi)
    {
        const fvPatchScalarField& pp = this->p_.boundaryField()[patchi];
        const fvPatchScalarField& pTt = this->Tt_.boundaryField()[patchi];
        
        fvPatchScalarField& phMix = composition().e().boundaryField()[patchi];
        fvPatchScalarField& phtMix = this->het_.boundaryField()[patchi];
        fvPatchScalarField& phvMix = this->hevMix_.boundaryField()[patchi];
        fvPatchScalarField& phelMix = this->heelMix_.boundaryField()[patchi];
        fvPatchScalarField& phvelMix = this->hevel().boundaryField()[patchi];
        fvPatchScalarField& pTvMix = this->Tv_.boundaryField()[patchi];
        fvPatchScalarField& pzetavMix = this->zetav_.boundaryField()[patchi];
        fvPatchScalarField ptotZetavMix = pzetavMix;
        
        //- Initialisations for patches
        phvMix = 0.0;
        phelMix = 0.0;
        if(not downgradeSingleTv) 
        {
            phvelMix = 0.0;
            pTvMix = 0.0;
        }
        pzetavMix = 0.0;
        ptotZetavMix = 0.0;
        
        if (pTt.fixesValue())
        // the temperature is fixed ... the energy is calculated at patches
        {
            phtMix = 0.0;
            
            forAll(Y, speciei)
            {
                const fvPatchScalarField& pY = Y[speciei].boundaryField()[patchi];
                const fvPatchScalarField& pX = X[speciei].boundaryField()[patchi];
                fvPatchScalarField& pTv = Tv[speciei].boundaryField()[patchi];
                fvPatchScalarField& phv = hv[speciei].boundaryField()[patchi];
                fvPatchScalarField& phel = hel[speciei].boundaryField()[patchi];
                fvPatchScalarField& phvel = hvel[speciei].boundaryField()[patchi];
                fvPatchScalarField& pzetav = zetav[speciei].boundaryField()[patchi];
                
                forAll(pTt, facei)
                {  
                    phtMix[facei] += pY[facei]*composition().HEt(speciei, pp[facei], pTt[facei]);
                    
                    if(pTv[facei] != 0.0)
                    {
                        if(downgradeSingleTemperature)
                        {
                            phv[facei] = composition().HEv(speciei, pp[facei], pTt[facei]);
                            phel[facei] = composition().HEel(speciei, pp[facei], pTt[facei]);
                            phvel[facei] = composition().HEvel(speciei, pp[facei], pTt[facei]);
                            phvelMix[facei] += pY[facei]*phvel[facei];
                            pzetav[facei] = composition().zetav(speciei, pp[facei], pTv[facei]);
                            pTvMix[facei] += pzetav[facei]*pTv[facei];
                        }
                        else if(downgradeSingleTv)
                        {
                            pTv[facei] = pTvMix[facei]; // NEW VINCENT 15/08/2016
                            phv[facei] = composition().HEv(speciei, pp[facei], pTt[facei]);
                            phel[facei] = composition().HEel(speciei, pp[facei], pTt[facei]);
                            phvel[facei] = composition().HEvel(speciei, pp[facei], pTv[facei]);
                            pzetav[facei] = composition().zetav(speciei, pp[facei], pTv[facei]);
                        }
                        else
                        {
                            phv[facei] = composition().HEv(speciei, pp[facei], pTt[facei]);
                            phel[facei] = composition().HEel(speciei, pp[facei], pTt[facei]);
                            phvel[facei] = composition().HEvel(speciei, pp[facei], pTv[facei]);
                            phvelMix[facei] += pY[facei]*phvel[facei];
                            pzetav[facei] = composition().zetav(speciei, pp[facei], pTv[facei]);
                            pTvMix[facei] += pzetav[facei]*pTv[facei];
                        }
                        
                        phvMix[facei] += pY[facei]*phv[facei];
                        phelMix[facei] += pY[facei]*phel[facei];
                        pzetavMix[facei] += pX[facei]*pzetav[facei];
                        ptotZetavMix[facei] += pzetav[facei];
                    }
                }//end faces loop
            }//end species loop
            
            phMix = phtMix + phvelMix;
            
            forAll(pTt, facei)
            {
                if(ptotZetavMix[facei] != 0.0)
                {
                    pTvMix[facei] /= ptotZetavMix[facei]; // NEW VINCENT 29/07/16  
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
                    
                    fvPatchScalarField& phvel = hvel[speciei].boundaryField()[patchi];
                    fvPatchScalarField& pTv = Tv[speciei].boundaryField()[patchi];
                    fvPatchScalarField& pzetav = zetav[speciei].boundaryField()[patchi];
                    
                    if (pTv[facei] != 0.0)
                    {
                        if(downgradeSingleTemperature)
                        {
                            pTv[facei] = pTt[facei];
                            phvel[facei] = composition().HEvel(speciei, pp[facei], pTv[facei]);
                            pzetav[facei] = composition().zetav(speciei, pp[facei], pTv[facei]);
                        }
                        else if(downgradeSingleTv)
                        {
                            phvel[facei] = composition().HEvel(speciei, pp[facei], pTv[facei]); // NEW VINCENT 27/01/2016 TODO
                            pzetav[facei] = composition().zetav(speciei, pp[facei], pTv[facei]);
                        }
                        else if (composition().vibTempAssociativity(speciei) == -1)
                        {
                            if(composition().noVibrationalTemp(speciei) == 1 or downgradeSingleVibMode)
                            {
                                if(pTt[facei] > vibrationalCutOffTemp)
                                {
                                    pTv[facei] = TvelEvels(speciei, phvel[facei], pp[facei], pTv[facei]);
                                    //composition().TvelHEvel(speciei, phvel[facei], pp[facei], pTv[facei]);
                                    pzetav[facei] = composition().zetav(speciei, pp[facei], pTv[facei]);
                                }
                                else
                                {
                                    pTv[facei] = pTt[facei];
                                    pzetav[facei] = composition().zetav(speciei, pp[facei], pTv[facei]);
                                }
                            }
                            else
                            {
                                /*pzetav[facei] = 0.0; TODO ONGOING WORK 
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
                        else // ions, electrons, and atoms if Eel is on
                        {
                            const fvPatchScalarField& pTvMol = Tv[composition().vibTempAssociativity(speciei)].boundaryField()[patchi];
                            pTv[facei] = pTvMol[facei];
                            pzetav[facei] = composition().zetav(speciei, pp[facei], pTv[facei]);
                        }
                    
                        phvMix[facei] += pY[facei]*phv[facei];
                        phelMix[facei] += pY[facei]*phel[facei];
                        phvelMix[facei] += pY[facei]*phvel[facei];
                        pTvMix[facei] += pzetav[facei]*pTv[facei];
                        pzetavMix[facei] += pX[facei]*pzetav[facei];
                        ptotZetavMix[facei] += pzetav[facei];
                    }
                }//end species loop
                
                if(ptotZetavMix[facei] != 0.0)
                {
                    pTvMix[facei] /= ptotZetavMix[facei];
                }
            }//end faces loop
        }
        
        fvPatchScalarField& ppsiMix = this->psi_.boundaryField()[patchi];
        ppsiMix = 0;
        
        forAll(Y, speciei)
        {
            const fvPatchScalarField& pX = X[speciei].boundaryField()[patchi];
            const fvPatchScalarField& pTv = Tv[speciei].boundaryField()[patchi];
            
            forAll(pTt, facei) 
            {
                if(composition().particleType(speciei) > 0) 
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
    
    scalarField& TtCells = this->Tt_.internalField();
    scalarField& TvCellsMix = this->Tv_.internalField(); // NEW VINCENT 15/08/2016
    scalarField& hvelCellsMix = this->hevel().internalField();
    scalarField& psiCellsMix = this->psi_.internalField();
       
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

        if(downgradeSingleTemperature)
        {
            TtCells[celli] = TEs
            (
                hCellsMix[celli],
                pCells[celli],
                TtCells[celli],
                YList
            );
        }
        else if(downgradeSingleTv)
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
        const volScalarField::GeometricBoundaryField wallPatches = this->Tt_.boundaryField(); // NEW VINCENT 23/08/2016
        const fvPatchScalarField& pp = this->p_.boundaryField()[patchi];
        const fvPatchScalarField& phMix = composition().e().boundaryField()[patchi];
        const fvPatchScalarField& phvelMix = this->hevel().boundaryField()[patchi]; // NEW VINCENT 16/08/2016
        fvPatchScalarField& phtMix = this->het_.boundaryField()[patchi];
        fvPatchScalarField& pTt = this->Tt_.boundaryField()[patchi];
        fvPatchScalarField& pTvMix = this->Tv_.boundaryField()[patchi]; // NEW VINCENT 16/08/2016
        
        if(pTt.fixesValue())
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
            
                if(downgradeSingleTemperature)
                {
                    pTt[facei] = TEs
                    (
                        phMix[facei],
                        pp[facei],
                        pTt[facei],
                        YList
                    );
                }
                else if(downgradeSingleTv) // NEW VINCENT 16/08/2016
                {
                    //Info << "phelo 1 " << endl;
                    pTt[facei] = TtEts
                    (
                        phtMix[facei],
                        pp[facei],
                        pTt[facei],
                        YList
                    );
                    
                    //Info << "phelo B " << patchi << tab << facei << endl;
                    pTvMix[facei] = TvelEvels
                    (
                        phvelMix[facei],
                        pp[facei],
                        pTvMix[facei],
                        YList
                    );
                    //Info << "phelo A " << patchi << tab << facei << endl;
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
        
        if(isA<wallFvPatch>(wallPatches[patchi].patch())) // NEW VINCENT 22/08/2016
        // to increase the stability of the 2T model in near-wall high-Kn number regions
        {
            forAll(pTt, facei)
            {
                if(pTt[facei] > ThighPatches_) pTt[facei] = ThighPatches_;
                else if(pTt[facei] < TlowPatches_) pTt[facei] = TlowPatches_;
            }
        }
    }//end patches loop
    
    
    //- Re-set of mixture quantities internalField
    if(not downgradeSingleTv) hvelCellsMix = 0.0;
    psiCellsMix = 0.0;

    //- Cells values calculation
    forAll(Y, speciei)
    {
        const scalarField& YCells = Y[speciei].internalField();
        const scalarField& XCells = X[speciei].internalField();
        scalarField& TvCells = Tv[speciei].internalField();
        scalarField& hvelCells = hvel[speciei].internalField();
            
        forAll(YCells, celli)
        {  
            if(TvCells[celli] != 0.0)
            {
                if(downgradeSingleTemperature)
                {
                    TvCells[celli] = TtCells[celli];
                    hvelCells[celli] = composition().HEvel(speciei, pCells[celli], TvCells[celli]);
                    hvelCellsMix[celli] += YCells[celli]*hvelCells[celli];
                }
                else if(downgradeSingleTv)
                {
                    TvCells[celli] = TvCellsMix[celli]; // NEW VINCENT 15/08/2016
                    hvelCells[celli] = composition().HEvel(speciei, pCells[celli], TvCells[celli]);
                }
                else if(composition().vibTempAssociativity(speciei) == -1) // NEW VINCENT 11/08/2016
                {
                    if(TtCells[celli] > vibrationalCutOffTemp)
                    {
                        if(YCells[celli] > miniYforSolvingEvEqn) // NEW VINCENT 08/02/2017 -> < means that evEqn not solved because Y_m too small
                        {
                            if(sign(hvelCells[celli]) == -1)
                            {
                                hvelCells[celli] = composition().HEvel(speciei, pCells[celli], vibrationalCutOffTemp); // NEW VINCENT 09/02/2017
                            }
                            //Info << composition().species()[speciei] << ": " << YCells[celli] << tab << hvelCells[celli];
                            TvCells[celli] = TvelEvels(speciei, hvelCells[celli], pCells[celli], TvCells[celli]);
                            //Info << tab << TvCells[celli] << endl;
                            //composition().TvelHEvel(speciei, hvelCells[celli], pCells[celli], TvCells[celli]);
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
            
            if(composition().particleType(speciei) > 0) 
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
        fvPatchScalarField& phtMix = this->het_.boundaryField()[patchi];
        fvPatchScalarField& phMix = composition().e().boundaryField()[patchi];
        fvPatchScalarField& phvelMix = this->hevel().boundaryField()[patchi];
        fvPatchScalarField& ppsiMix = this->psi_.boundaryField()[patchi];
        
        //- Re-set of mixture quantities boundaryField
        forAll(pTt, facei)
        {
            if(not downgradeSingleTv) phvelMix[facei] = 0.0;
            ppsiMix[facei] = 0.0;
        }
        
        if(pTt.fixesValue())
        // the temperature is fixed ... the energy is calculated at patches
        {
            forAll(Y, speciei)
            {
                const fvPatchScalarField& pY = Y[speciei].boundaryField()[patchi];
                fvPatchScalarField& pTv = Tv[speciei].boundaryField()[patchi];
                fvPatchScalarField& phvel = hvel[speciei].boundaryField()[patchi];
                
                forAll(pTt, facei)
                {  
                    if(pTv[facei] != 0.0)
                    {
                        if(downgradeSingleTemperature)
                        {
                            phvel[facei] = composition().HEvel(speciei, pp[facei], pTt[facei]);
                            phvelMix[facei] += pY[facei]*phvel[facei];
                        }
                        else if(downgradeSingleTv)
                        {
                            pTv[facei] = pTvMix[facei]; // NEW VINCENT 15/08/2016
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
                fvPatchScalarField& phvel = hvel[speciei].boundaryField()[patchi];
                fvPatchScalarField& pTv = Tv[speciei].boundaryField()[patchi];
                
                forAll(pTt, facei)
                {
                    if (pTv[facei] != 0.0)
                    {
                        if(downgradeSingleTemperature)
                        {
                            pTv[facei] = pTt[facei];
                            phvel[facei] = composition().HEvel(speciei, pp[facei], pTv[facei]);
                            phvelMix[facei] += pY[facei]*phvel[facei];
                        }
                        else if(downgradeSingleTv)
                        {
                            pTv[facei] = pTvMix[facei]; // NEW VINCENT 15/08/2016
                            phvel[facei] = composition().HEvel(speciei, pp[facei], pTv[facei]);
                        }
                        else if (composition().vibTempAssociativity(speciei) == -1) // NEW VINCENT 11/08/2016
                        {
                            if(pTt[facei] > vibrationalCutOffTemp)
                            {
                                if(pY[facei] > miniYforSolvingEvEqn) // NEW VINCENT 08/02/2017 -> < means that evEqn not solved because Y_m too small
                                {
                                    if(sign(phvel[facei]) == -1) // NEW VINCENT 09/02/2017
                                    {
                                        phvel[facei] = composition().HEvel(speciei, pp[facei], vibrationalCutOffTemp);
                                    }
                                    //Info << composition().species()[speciei] << ": " << pY[facei] << tab << phvel[facei];
                                    pTv[facei] = TvelEvels(speciei, phvel[facei], pp[facei], pTv[facei]);
                                    //Info << tab << pTv[facei] << endl;
                                    //composition().TvelHEvel(speciei, phvel[facei], pp[facei], pTv[facei]);
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
                if(composition().particleType(speciei) > 0) 
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
    
    //- Declarations
    const bool downgradeSingleTv = this->downgradeSingleTv();
    
    const scalarField& pCells = this->p_.internalField();
    const scalarField& TvCellsMix = this->Tv_.internalField();
    PtrList<Foam::volScalarField>& Tv = composition().Tv();
    PtrList<Foam::volScalarField>& hvel = composition().hevel();
       
    //- Cells values
    forAll(Tv, speciei)
    {  
        scalarField& TvCells = Tv[speciei].internalField();
        scalarField& hvelCells = hvel[speciei].internalField();

        if(TvCells[0] != 0.0)
        {
            if(downgradeSingleTv)
            {
                forAll(pCells, celli)
                {
                    TvCells[celli] = TvCellsMix[celli];
                    hvelCells[celli] = composition().HEvel(speciei, pCells[celli], TvCells[celli]);
                }
            }
            else if(composition().vibTempAssociativity(speciei) != -1) // NEW VINCENT 04/08/2016
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
            fvPatchScalarField& pTv = Tv[speciei].boundaryField()[patchi];
            const fvPatchScalarField& pTvMix = this->Tv_.boundaryField()[patchi];
            
            if(pTv.size() != 0) // NEW VINCENT 04/08/2016
            {    
                if(pTv[0] != 0.0)
                {
                    const fvPatchScalarField& pp = this->p_.boundaryField()[patchi];
                    fvPatchScalarField& phvel = hvel[speciei].boundaryField()[patchi];
                        
                    if(downgradeSingleTv)
                    {
                        forAll(pTv, facei)
                        {
                            pTv[facei] = pTvMix[facei];
                            phvel[facei] = composition().HEvel(speciei, pp[facei], pTv[facei]);
                        }
                    }
                    else if(composition().vibTempAssociativity(speciei) != -1) // NEW VINCENT 04/08/2016
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


void Foam::rho2ReactionThermo::limitVelocityAtWallBoundary(volVectorField::GeometricBoundaryField& Ubdry)
{
    const volScalarField::GeometricBoundaryField wallPatches = this->Tt_.boundaryField();
    
    //- Patch values
    forAll(Ubdry, patchi)
    {
        if(isA<wallFvPatch>(wallPatches[patchi].patch()))
        // to increase the stability of the 2T model in near-wall high Kn-number regions
        {
            fvPatchVectorField& pU = Ubdry[patchi];
            
            forAll(pU, facei)
            {
                if(pU[facei].component(0) > UhighPatches_) pU[facei].component(0) = UhighPatches_;
                if(pU[facei].component(1) > UhighPatches_) pU[facei].component(1) = UhighPatches_;
                if(pU[facei].component(2) > UhighPatches_) pU[facei].component(2) = UhighPatches_;
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::rho2ReactionThermo::rho2ReactionThermo
(
    const fvMesh& mesh,
    const word& phaseName
)
:
    rho2Thermo(mesh, phaseName),
    partialThermoName_(transportToTypedef(word(subDict("thermoType").lookup("transport"))))
{
    if(this->downgradeSingleTv_)
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
    
    if(isDict("temperatureBounds"))
    {
        Tlow_ = subDict("temperatureBounds").lookupOrDefault<scalar>("Tlow", 100);
        Thigh_ = subDict("temperatureBounds").lookupOrDefault<scalar>("Thigh", 40000);
        TlowPatches_ = subDict("temperatureBounds").lookupOrDefault<scalar>("TlowPatches", Tlow_);
        ThighPatches_ = subDict("temperatureBounds").lookupOrDefault<scalar>("ThighPatches", Thigh_);
    }
    else
    {
        Tlow_ = 100;
        Thigh_ = 40000;
        TlowPatches_ = Tlow_;
        ThighPatches_ = Thigh_;
    }
    
    if(isDict("velocityBounds"))
    {
        UhighPatches_ = subDict("velocityBounds").lookupOrDefault<scalar>("UhighPatches", Foam::GREAT);
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


// NEW VINCENT ****************************************************************
// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::rho2ReactionThermo& Foam::rho2ReactionThermo::lookup2ReactionThermo
(
    const fvPatchScalarField& pf
)
{
    return pf.db().lookupObject<rho2ReactionThermo>(dictName); // NOTE VINCENT 20/02/2016: adapted from basicThermo.C
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
    PtrList<volScalarField>& fixedY,
    volScalarField& fixedTt,
    PtrList<volScalarField>& fixedTv,
    volScalarField& fixedRho,
    word& regionName
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
// END NEW VINCENT ************************************************************


// ************************************************************************* //
