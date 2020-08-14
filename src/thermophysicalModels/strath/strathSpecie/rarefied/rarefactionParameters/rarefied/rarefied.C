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

#include "rarefied.H"
#include "fvc.H"

// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

template<class ThermoType>
void Foam::rarefied<ThermoType>::updateCoefficients()
{
    mfpModel_().update();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ThermoType>
Foam::rarefied<ThermoType>::rarefied
(
    rho2ReactionThermo& thermo,
    const compressibleTurbulenceModel& turbulence
)
:
    rarefactionParameter(thermo, turbulence),

    speciesThermo_
    (
        dynamic_cast<const multi2ComponentMixture<ThermoType>&>
            (this->thermo_).speciesData()
    ),

    miniXs_(1e-4)
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class ThermoType>
void Foam::rarefied<ThermoType>::correct(const volVectorField& U)
{
    if(computeRarefaction_)
    {
        updateCoefficients();

        const volScalarField& T = thermo_.T();
        const volScalarField& Tv = thermo_.Tv();
        const volScalarField magU (mag(U));
        volScalarField c(sqrt(thermo_.Cp_t()/thermo_.Cv_t()/thermo_.psi()));

        const volScalarField& p = thermo_.p();
        const volScalarField& rho = thermo_.rho();
        const volScalarField& muMix = thermo_.mu();

        const scalarField& pCells = p.internalField();
        const scalarField& TCells = T.internalField();
        const scalarField& TvCells = Tv.internalField();
        const scalarField& UCells = magU.internalField();
        const scalarField& cCells = c.internalField();
        const scalarField& rhoCells = rho.internalField();
        const scalarField& muMixCells = muMix.internalField();

        scalarField& mfpMixCells = mfpMix_.primitiveFieldRef();
        volScalarField innerSum = mfpMix_, outerSum = mfpMix_; // NEW VINCENT 03/08/2016
        volScalarField nDmix = thermo_.composition().nD()[0]; // NEW VINCENT 03/08/2016

        volScalarField tMfpMix = mfpMix_;

        //- Initialisation
        forAll(rhoCells, celli)
        {
            innerSum[celli] = 0.0;
            outerSum[celli] = 0.0;
            nDmix[celli] = 0.0;
            tMfpMix[celli] = 0.0;
        }

        forAll(rho.boundaryField(), patchi)
        {
            forAll(rho.boundaryField()[patchi], facei)
            {
                innerSum.boundaryFieldRef()[patchi][facei] = 0.0;
                outerSum.boundaryFieldRef()[patchi][facei] = 0.0;
                nDmix.boundaryFieldRef()[patchi][facei] = 0.0;
                tMfpMix.boundaryFieldRef()[patchi][facei] = 0.0;
            }
        }

        forAll(species(), speciei)
        {
            const volScalarField& X = thermo_.composition().X()[speciei];
            const volScalarField& pD = thermo_.composition().pD()[speciei];
            volScalarField& mfp = mfpModel_().mfp(speciei); // NOTE VINCENT: non-constant declaration because it has to be multiplied by mu_i/rho_i

            const scalarField& XCells = X.internalField();
            const scalarField& pDCells = pD.internalField();

            scalarField& mfpCells = mfp.primitiveFieldRef();

            if(oldMfpDefinition_) // NEW VINCENT 03/08/2016
            {
                tMfpMix += X*mfp;
            }

            forAll(rhoCells, celli)
            {
                if(XCells[celli] > miniXs_)
                {
                    mfpCells[celli] *= mu(speciei, pCells[celli], TCells[celli])/pDCells[celli];
                }
                else
                {
                    mfpCells[celli] = Foam::GREAT;
                }
            }

            forAll(rho.boundaryField(), patchi)
            {
                const fvPatchScalarField& pX = X.boundaryField()[patchi];
                const fvPatchScalarField& ppD = pD.boundaryField()[patchi];
                const fvPatchScalarField& pp = p.boundaryField()[patchi];
                const fvPatchScalarField& pT = T.boundaryField()[patchi];
                fvPatchScalarField& pmfp = mfp.boundaryFieldRef()[patchi];

                forAll(pp, facei)
                {
                    if(pX[facei] > miniXs_)
                    {
                        pmfp[facei] *= mu(speciei, pp[facei], pT[facei])/ppD[facei];
                    }
                    else
                    {
                        pmfp[facei] = Foam::GREAT;
                    }
                }
            }


            if(not oldMfpDefinition_ /*and species().size() > 1*/) // NEW VINCENT 03/08/2016
            {
                const volScalarField& nDi = thermo_.composition().nD()[speciei];
                const scalarField& nDiCells = nDi.internalField();

                //- Initialisation
                forAll(rhoCells, celli)
                {
                    innerSum[celli] = 0.0;
                }

                forAll(rho.boundaryField(), patchi)
                {
                    forAll(rho.boundaryField()[patchi], facei)
                    {
                        innerSum.boundaryFieldRef()[patchi][facei] = 0.0;
                    }
                }

                //- Calculation
                forAll(species(), speciej)
                {
                    const volScalarField& nDj = thermo_.composition().nD()[speciej];
                    const scalarField& nDjCells = nDj.internalField();

                    forAll(T, celli)
                    {
                        innerSum[celli] += max(nDjCells[celli]*innerQuantity(speciei, speciej, TCells[celli]), Foam::VSMALL);
                    }

                    forAll(T.boundaryField(), patchi)
                    {
                        const fvPatchScalarField& pT = T.boundaryField()[patchi];
                        const fvPatchScalarField& pnDj = nDj.boundaryField()[patchi];
                        fvPatchScalarField& pInnerSum = innerSum.boundaryFieldRef()[patchi];

                        forAll(pT, facei)
                        {
                            pInnerSum[facei] += max(pnDj[facei]*innerQuantity(speciei, speciej, pT[facei]), Foam::VSMALL);
                        }
                    }
                }// end speciej loop

                forAll(T, celli)
                {
                    outerSum[celli] += nDiCells[celli]/innerSum[celli];
                }

                forAll(nDi.boundaryField(), patchi)
                {
                    const fvPatchScalarField& pnDi = nDi.boundaryField()[patchi];
                    const fvPatchScalarField& pInnerSum = innerSum.boundaryField()[patchi];
                    fvPatchScalarField& pOuterSum = outerSum.boundaryFieldRef()[patchi];

                    forAll(pnDi, facei)
                    {
                        pOuterSum[facei] += pnDi[facei]/pInnerSum[facei];
                    }
                }

                nDmix += nDi;
            }
        }

        const scalarField& tMfpMixCells = tMfpMix.internalField();
        scalarField& KnovCells = Knov_.primitiveFieldRef();
        scalarField& KnGLLCells = KnGLL_.primitiveFieldRef();
        scalarField& KnGLLCellsRho = KnsGLL_[0].primitiveFieldRef();
        scalarField& KnGLLCellsT = KnsGLL_[1].primitiveFieldRef();
        scalarField& KnGLLCellsU = KnsGLL_[2].primitiveFieldRef();

        const volScalarField gradRho = mag(fvc::grad(thermo_.rho()));
        const scalarField& gradRhoCells = gradRho.internalField();
        const volScalarField gradT = mag(fvc::grad(thermo_.T()));
        const volScalarField gradTv = mag(fvc::grad(thermo_.Tv()));
        const scalarField& gradTCells = gradT.internalField();
        const scalarField& gradTvCells = gradTv.internalField();
        const volScalarField gradU = mag(fvc::grad(magU));
        const scalarField& gradUCells = gradU.internalField();

        forAll(rhoCells, celli)
        {
            if(oldMfpDefinition_) // NEW VINCENT 03/08/2016
            {
                mfpMixCells[celli] = tMfpMixCells[celli]*muMixCells[celli]/rhoCells[celli];
            }
            else /*if(species().size() > 1)*/
            {
                mfpMixCells[celli] = outerSum[celli]/nDmix[celli];
            }
            /*else
            {
                mfpMixCells[celli] = mfpModel_().mfp(0)[celli];
            }*/

            KnGLLCellsRho[celli] = mfpMixCells[celli]/rhoCells[celli]*gradRhoCells[celli];
            KnGLLCellsU[celli] = mfpMixCells[celli]/max(max(cCells[celli], UCells[celli]), Foam::SMALL)*gradUCells[celli]; // NOTE VINCENT: to be able to run heat baths
            KnovCells[celli] = mfpMixCells[celli]/characteristicLength_;
        }

        if(TvCells[0] != 0.0)
        {
            forAll(rhoCells, celli)
            {
                KnGLLCellsT[celli] = mfpMixCells[celli]*max(gradTCells[celli]/TCells[celli], gradTvCells[celli]/TvCells[celli]);
                KnGLLCells[celli] = max(max(KnGLLCellsRho[celli], KnGLLCellsT[celli]), KnGLLCellsU[celli]);
            }
        }
        else
        {
            forAll(rhoCells, celli)
            {
                KnGLLCellsT[celli] = mfpMixCells[celli]*gradTCells[celli]/TCells[celli];
                KnGLLCells[celli] = max(max(KnGLLCellsRho[celli], KnGLLCellsT[celli]), KnGLLCellsU[celli]);
            }
        }

        forAll(rho.boundaryField(), patchi)
        {
            const fvPatchScalarField& pmuMix = muMix.boundaryField()[patchi];
            const fvPatchScalarField& ptMfpMix = tMfpMix.boundaryField()[patchi];
            const fvPatchScalarField& prho = rho.boundaryField()[patchi];
            const fvPatchScalarField& pT = T.boundaryField()[patchi];
            const fvPatchScalarField& pTv = Tv.boundaryField()[patchi];
            const fvPatchScalarField& pU = magU.boundaryField()[patchi];
            const fvPatchScalarField& pc = c.boundaryField()[patchi];
            const fvPatchScalarField& pgradRho = gradRho.boundaryField()[patchi];
            const fvPatchScalarField& pgradT = gradT.boundaryField()[patchi];
            const fvPatchScalarField& pgradTv = gradTv.boundaryField()[patchi];
            const fvPatchScalarField& pgradU = gradU.boundaryField()[patchi];

            const fvPatchScalarField& pOuterSum = outerSum.boundaryField()[patchi]; // NEW VINCENT 03/08/2016
            const fvPatchScalarField& pnDmix = nDmix.boundaryField()[patchi]; // NEW VINCENT 03/08/2016

            fvPatchScalarField& pmfpMix = mfpMix_.boundaryFieldRef()[patchi];
            fvPatchScalarField& pKnov = Knov_.boundaryFieldRef()[patchi];
            fvPatchScalarField& pKnGLL = KnGLL_.boundaryFieldRef()[patchi];
            fvPatchScalarField& pKnGLLRho = KnsGLL_[0].boundaryFieldRef()[patchi];
            fvPatchScalarField& pKnGLLT = KnsGLL_[1].boundaryFieldRef()[patchi];
            fvPatchScalarField& pKnGLLU = KnsGLL_[2].boundaryFieldRef()[patchi];

            forAll(prho, facei)
            {
                if(oldMfpDefinition_) // NEW VINCENT 03/08/2016
                {
                    pmfpMix[facei] = ptMfpMix[facei]*pmuMix[facei]/prho[facei];
                }
                else /*if(species().size() > 1)*/
                {
                    pmfpMix[facei] = pOuterSum[facei]/pnDmix[facei];
                }
                /*else
                {
                    pmfpMix[facei] = mfpModel_().mfp(0).boundaryField()[patchi][facei];
                }*/

                pKnGLLRho[facei] = pmfpMix[facei]/prho[facei]*pgradRho[facei];
                if(pTv[0] != 0.0)
                {
                    pKnGLLT[facei] = pmfpMix[facei]*max(pgradT[facei]/pT[facei], pgradTv[facei]/pTv[facei]);
                }
                else
                {
                    pKnGLLT[facei] = pmfpMix[facei]*pgradT[facei]/pT[facei];
                }
                pKnGLLU[facei] = pmfpMix[facei]/max(max(pU[facei], pc[facei]), Foam::SMALL)*pgradU[facei];
                pKnGLL[facei] = max(max(pKnGLLRho[facei], pKnGLLT[facei]), pKnGLLU[facei]);
                pKnov[facei] = pmfpMix[facei]/characteristicLength_;
            }
        }
    }
    else if(computeMfpBoundaries_)
    {
        if(oldMfpDefinition_)
        {
            updateCoefficients();

            const volScalarField& T = thermo_.T();
            const volScalarField& p = thermo_.p();
            const volScalarField& rho = thermo_.rho();
            const volScalarField& muMix = thermo_.mu();

            volScalarField tMfpMix = mfpMix_*0.0;

            forAll(species(), speciei)
            {
                const volScalarField& X = thermo_.composition().X()[speciei];
                const volScalarField& pD = thermo_.composition().pD()[speciei];
                volScalarField& mfp = mfpModel_().mfp(speciei); // NOTE VINCENT: non-constant declaration because it has to be multiplied by mu_i/rho_i

                tMfpMix.boundaryFieldRef() += X.boundaryField()*mfp.boundaryField();

                forAll(rho.boundaryField(), patchi)
                {
                    const fvPatchScalarField& pX = X.boundaryField()[patchi];
                    const fvPatchScalarField& ppD = pD.boundaryField()[patchi];
                    const fvPatchScalarField& pp = p.boundaryField()[patchi];
                    const fvPatchScalarField& pT = T.boundaryField()[patchi];
                    fvPatchScalarField& pmfp = mfp.boundaryFieldRef()[patchi];

                    forAll(pp, facei)
                    {
                        if(pX[facei] > miniXs_)
                        {
                            pmfp[facei] *= mu(speciei, pp[facei], pT[facei])/ppD[facei];
                        }
                        else
                        {
                            pmfp[facei] = Foam::GREAT;
                        }
                    }
                }
            }

            forAll(rho.boundaryField(), patchi)
            {
                const fvPatchScalarField& pmuMix = muMix.boundaryField()[patchi];
                const fvPatchScalarField& ptMfpMix = tMfpMix.boundaryField()[patchi];
                const fvPatchScalarField& prho = rho.boundaryField()[patchi];

                fvPatchScalarField& pmfpMix = mfpMix_.boundaryFieldRef()[patchi];

                forAll(prho, facei)
                {
                    pmfpMix[facei] = ptMfpMix[facei]*pmuMix[facei]/prho[facei];
                }
            }
        }
        else
        {
            const volScalarField& Tt = thermo_.T();
            const scalarField& TtCells = Tt.internalField();

            volScalarField innerSum = mfpMix_, outerSum = mfpMix_;
            volScalarField nDmix = thermo_.composition().nD()[0];

            //- Initialisation
            forAll(TtCells, celli)
            {
                innerSum[celli] = 0.0;
                outerSum[celli] = 0.0;
                nDmix[celli] = 0.0;
            }

            forAll(Tt.boundaryField(), patchi)
            {
                forAll(Tt.boundaryField()[patchi], facei)
                {
                    innerSum.boundaryFieldRef()[patchi][facei] = 0.0;
                    outerSum.boundaryFieldRef()[patchi][facei] = 0.0;
                    nDmix.boundaryFieldRef()[patchi][facei] = 0.0;
                }
            }

            forAll(species(), speciei)
            {
                const volScalarField& nDi = thermo_.composition().nD()[speciei];

                //- Initialisation
                forAll(TtCells, celli)
                {
                    innerSum[celli] = 0.0;
                }

                forAll(Tt.boundaryField(), patchi)
                {
                    forAll(Tt.boundaryField()[patchi], facei)
                    {
                        innerSum.boundaryFieldRef()[patchi][facei] = 0.0;
                    }
                }

                //- Calculation
                forAll(species(), speciej)
                {
                    const volScalarField& nDj = thermo_.composition().nD()[speciej];

                    forAll(Tt.boundaryField(), patchi)
                    {
                        const fvPatchScalarField& pTt = Tt.boundaryField()[patchi];
                        const fvPatchScalarField& pnDj = nDj.boundaryField()[patchi];
                        fvPatchScalarField& pInnerSum = innerSum.boundaryFieldRef()[patchi];

                        forAll(pTt, facei)
                        {
                            pInnerSum[facei] += max(pnDj[facei]*innerQuantity(speciei, speciej, pTt[facei]), Foam::VSMALL);
                        }
                    }
                }

                forAll(nDi.boundaryField(), patchi)
                {
                    const fvPatchScalarField& pnDi = nDi.boundaryField()[patchi];
                    const fvPatchScalarField& pInnerSum = innerSum.boundaryField()[patchi];
                    fvPatchScalarField& pOuterSum = outerSum.boundaryFieldRef()[patchi];

                    forAll(pnDi, facei)
                    {
                        pOuterSum[facei] += pnDi[facei]/pInnerSum[facei];
                    }
                }

                nDmix.boundaryFieldRef() += nDi.boundaryField();
            }

            forAll(mfpMix_.boundaryField(), patchi)
            {
                const fvPatchScalarField& pOuterSum = outerSum.boundaryField()[patchi];
                const fvPatchScalarField& pnDmix = nDmix.boundaryField()[patchi];
                fvPatchScalarField& pmfpMix = mfpMix_.boundaryFieldRef()[patchi];

                forAll(pmfpMix, facei)
                {
                    pmfpMix[facei] = pOuterSum[facei]/pnDmix[facei];
                }
            }
        }
    }
}


template<class ThermoType>
void Foam::rarefied<ThermoType>::write()
{
    if(computeRarefaction_)
    {
        if (writeMfpSpecies_)
        {
            forAll(species(), speciei)
            {
                mfpModel_().mfp(speciei).write();
            }
        }

        if (writeMfpMixture_)
        {
            mfpMix_.write();
        }

        if (writeKnGLL_)
        {
            KnGLL_.write();
        }

        if (writeKnGLLComponents_)
        {
            forAll(KnsGLL_, i)
            {
                KnsGLL_[i].write();
            }
        }

        if (writeKnOv_)
        {
            Knov_.write();
        }
    }
}


template<class ThermoType>
bool Foam::rarefied<ThermoType>::read()
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
